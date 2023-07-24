version 1.0

workflow Soma {
    input {
        File InputSpreadSheet
        File SampleSheet
        # sample sheet has this structure:
        # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE [R1] [R2]

        File? DemuxSampleSheet
        String? IlluminaDir
        String? DragenEnv

        Boolean DataTransfer
        Boolean RmRunDir

        String XferLabel
        String OutputDir
        String Queue
        String JobGroup
        String DragenQueue
        String DragenDockerImage

        String SomaRepo
        String CoverageBed  = SomaRepo + "/accessory_files/SOMA.all.bed"
        String CovLevels = "100,500,1000,1500"
        String HaplotectBed = SomaRepo + "/accessory_files/SOMA.haplotect.bed"
        String QC_py        = SomaRepo + "/scripts/QC_metrics.py"
    }

    String DragenReference = "/storage1/fs1/gtac-mgi/Active/CLE/reference/dragen_hg38"
    String Reference       = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.fa"
    String ReferenceDict   = "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/refdata/hg38/all_sequences.dict"

    String DemuxFastqDir = "/storage1/fs1/gtac-mgi/Active/CLE/assay/SOMA/demux_fastq"

    Int readfamilysize  = 1

    if (defined(DemuxSampleSheet)){
        call dragen_demux {
            input: Dir=IlluminaDir,
                   OutputDir=OutputDir,
                   SampleSheet=DemuxSampleSheet,
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call prepare_samples {
            input: SampleSheet=SampleSheet,
                   Fastq1=dragen_demux.read1,
                   Fastq2=dragen_demux.read2,
                   queue=Queue,
                   jobGroup=JobGroup
        }
    }

    Array[Array[String]] inputData = read_tsv(select_first([prepare_samples.sample_sheet,SampleSheet]))

    # the inputdata should be: index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE read1path read2path
    scatter (samples in inputData){
        call dragen_align {
            input: DragenRef=DragenReference,
                   fastq1=samples[7],
                   fastq2=samples[8],
                   Name=samples[1],
                   RG=samples[3] + '.' + samples[4] + '.' + samples[0],
                   SM=samples[6],
                   LB=samples[5] + '.' + samples[0],
                   readfamilysize=readfamilysize,
                   CoverageBed=CoverageBed,
                   CovLevels=CovLevels,
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call run_haplotect {
            input: refFasta=Reference,
                   refDict=ReferenceDict,
                   Cram=dragen_align.cram,
                   CramIndex=dragen_align.crai,
                   Bed=HaplotectBed,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call gather_files {
            input: OutputFiles=[run_haplotect.out_file,
                   run_haplotect.sites_file],
                   OutputDir=OutputDir,
                   SubDir=samples[1] + '_' + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
        }
    }

    call batch_qc {
        input: order_by=gather_files.done,
               InputSpreadSheet=InputSpreadSheet,
               BatchDir=OutputDir,
               QC_py=QC_py,
               queue=Queue,
               jobGroup=JobGroup
    }

    if (defined(DemuxSampleSheet)){
        call move_demux_fastq {
            input: order_by=gather_files.done,
                   Batch=basename(OutputDir),
                   DemuxFastqDir=DemuxFastqDir,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        if (DataTransfer) {
            call data_transfer {
                input: order_by=move_demux_fastq.done,
                QcFile=batch_qc.QC_file,
                XferLabel=XferLabel,
                BatchFastqDir= DemuxFastqDir + '/' + basename(OutputDir),
                queue=Queue,
                jobGroup=JobGroup
            }
        }

        if (RmRunDir) {
            call remove_rundir {
                input: order_by=move_demux_fastq.done,
                       rundir=IlluminaDir,
                       queue=DragenQueue,
                       jobGroup=JobGroup
            }
        }
    }
}


task dragen_demux {
     input {
         String? Dir
         String OutputDir
         File? SampleSheet
         String? DragenEnv
         String DragenDockerImage
         String jobGroup
         String queue
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/Soma/"
     String LocalFastqDir = StagingDir + "demux_fastq/" + batch
     String LocalReportDir = LocalFastqDir + "/Reports"
     String LocalSampleSheet = StagingDir + "sample_sheet/" + batch + '.csv'
     String log = StagingDir + "log/" + batch + "_demux.log"
     String DemuxReportDir = OutputDir + "/dragen_demux_reports"

     command <<<
         /bin/cp ~{SampleSheet} ~{LocalSampleSheet} && \
         /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ~{LocalSampleSheet} --bcl-input-directory ~{Dir} --output-directory ~{LocalFastqDir} &> ~{log} && \
         /bin/ls ~{LocalFastqDir}/*_R1_001.fastq.gz > Read1_list.txt && \
         /bin/ls ~{LocalFastqDir}/*_R2_001.fastq.gz > Read2_list.txt && \
         /bin/mv ~{log} ./ && \
         /bin/rm -f ~{LocalSampleSheet} && \
         /bin/cp -r ~{LocalReportDir} ~{DemuxReportDir}
     >>>

     runtime {
         docker_image: DragenDockerImage
         cpu: "20"
         memory: "200 G"
         dragen_env: DragenEnv
         queue: queue
         job_group: jobGroup
     }
     output {
         File read1 = "Read1_list.txt"
         File read2 = "Read2_list.txt"
     }
}

task prepare_samples {
     input {
         File SampleSheet
         String Fastq1
         String Fastq2
         String jobGroup
         String queue
     }

     command <<<
         /bin/cp ~{Fastq1} 1.tmp.txt
         /bin/cp ~{Fastq2} 2.tmp.txt
         /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
             chomp @r1; close R1;\
             open(R2,"2.tmp.txt"); @r2 = <R2>; \
             chomp @r2; close R2; \
             open(SS,"~{SampleSheet}");
             while(<SS>){
                 chomp;
                 my @l = split("\t",$_);
                 my $s = $l[1].'_';
                 my $r1 = (grep /$s/, @r1)[0];
                 my $r2 = (grep /$s/, @r2)[0];
                 print join("\t",@l,$r1,$r2),"\n";
             }
             close SS;' > sample_sheet.txt
     >>>

     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/genome/lims-compute-xenial:1)"
         cpu: "1"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }

     output {
         File sample_sheet = "sample_sheet.txt"
     }
}

task dragen_align {
     input {
         String Name
         String DragenRef
         String fastq1
         String fastq2
         String RG
         String SM
         String LB
         String CoverageBed
         String CovLevels
         String OutputDir
         String SubDir
         String DragenDockerImage
         String jobGroup
         String queue

         String? DragenEnv
         Int readfamilysize
     }

     String batch = basename(OutputDir)
     String StagingDir = "/staging/runs/Soma/"
     String LocalAlignDir = StagingDir + "align/" + batch
     String LocalSampleDir = LocalAlignDir + "/" + SubDir
     String log = StagingDir + "log/" + Name + "_align.log"

     String outdir = OutputDir + "/" + SubDir
     String dragen_outdir = outdir + "/dragen"


     command {
         if [ ! -d "${LocalAlignDir}" ]; then
             /bin/mkdir ${LocalAlignDir}
         fi

         /bin/mkdir ${LocalSampleDir} && \
         /bin/mkdir ${outdir} && \
         /opt/edico/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} --enable-map-align true --enable-sort true --enable-map-align-output true --enable-variant-caller=true --vc-enable-umi-solid true --vc-combine-phased-variants-distance 3 --vc-enable-orientation-bias-filter true --vc-enable-triallelic-filter false --vc-target-bed ${CoverageBed} --gc-metrics-enable=true --qc-coverage-ignore-overlaps=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 cov_report --qc-coverage-region-1-thresholds ${CovLevels} --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads ${readfamilysize} --umi-metrics-interval-file ${CoverageBed} --output-dir ${LocalSampleDir} --output-file-prefix ${Name} --output-format CRAM &> ${log} && \
         /bin/mv ${log} ./ && \
         /bin/mv ${LocalSampleDir} ${dragen_outdir}
     }

     runtime {
         docker_image: DragenDockerImage 
         cpu: "20"
         memory: "200 G"
         dragen_env: DragenEnv
         queue: queue
         job_group: jobGroup
     }

     output {
         File cram = "${dragen_outdir}/${Name}_tumor.cram"
         File crai = "${dragen_outdir}/${Name}_tumor.cram.crai"
     }
}

task run_haplotect {
     input {
         String Cram
         String CramIndex
         String Bed
         String Name
         String refDict
         String refFasta
         String queue
         String jobGroup

         Int? MinReads
     }
     command <<<
         /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ~{Bed} > /tmp/pos.bed && \
         /usr/local/openjdk-8/bin/java -Xmx6g \
         -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
         -I ~{Cram} -R ~{refFasta} --sequence-dictionary ~{refDict} \
         -mmq 20 -mbq 20 -max-depth-per-sample 10000 -gstol 0.001 -mr ~{default=10 MinReads} \
         -htp ~{Bed} -L /tmp/pos.bed -outPrefix ~{Name}
     >>>

     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/mgi-cle/haplotect:0.3)"
         cpu: "1"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File out_file = "${Name}.haplotect.txt"
         File sites_file = "${Name}.haplotectloci.txt"
     }
}

task gather_files {
     input {
         Array[String] OutputFiles
         String OutputDir
         String? SubDir
         String jobGroup
         String queue
     }
     command {
         if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
             mkdir ${OutputDir}/${SubDir}
         fi
         /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
     }
     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/genome/lims-compute-xenial:1)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task batch_qc {
     input {
         Array[String] order_by
         String InputSpreadSheet
         String BatchDir
         String QC_py
         String queue
         String jobGroup
     }
     String batch = basename(BatchDir)

     command {
         if [ -n "$(/bin/ls -d ${BatchDir}/H_*)" ]; then
             /bin/chmod -R 777 ${BatchDir}/H_*
         fi

         /usr/bin/python3 ${QC_py} ${InputSpreadSheet} ${BatchDir}
     }
     runtime {
         docker_image: "docker1(catgumag/pandas-scibioxl:20220107)"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File QC_file = "${BatchDir}/${batch}_Genoox.xlsx"
     }
}

task move_demux_fastq {
     input {
         Array[String] order_by
         String Batch
         String DemuxFastqDir
         String queue
         String jobGroup
     }

     String LocalDemuxFastqDir = "/staging/runs/Soma/demux_fastq/" + Batch

     command {
         if [ -d "${LocalDemuxFastqDir}" ]; then
             /bin/mv ${LocalDemuxFastqDir} ${DemuxFastqDir}
         fi
     }
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task remove_rundir {
     input {
         String order_by
         String? rundir
         String queue
         String jobGroup
     }
     command {
         if [ -d "${rundir}" ]; then
             /bin/rm -Rf ${rundir}
         fi
     }
     runtime {
         docker_image: "docker1(ubuntu:xenial)"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}

task data_transfer {
     input {
         String order_by
         String QcFile
         String XferLabel
         String BatchFastqDir
         String queue
         String jobGroup
     }
     
     command {
         set -eo pipefail && \
         /bin/mkdir xfer_staging && \
         /bin/cp ${QcFile} xfer_staging && \
         /bin/cp ${BatchFastqDir}/H_*.fastq.gz xfer_staging && \
         /usr/local/bin/aws s3 cp xfer_staging s3://genoox-upload-wustl/gtacmgi/${XferLabel} --recursive && \
         /usr/bin/touch done.txt && \
         /usr/local/bin/aws s3 cp done.txt s3://genoox-upload-wustl/gtacmgi/${XferLabel}
     }
     runtime {
         docker_image: "docker1(registry.gsc.wustl.edu/dataxfer/data-transfer-aws)"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
