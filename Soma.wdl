version 1.0

workflow Soma {
    input {
        File SampleSheet
        # sample sheet has this structure:
        # index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE [R1] [R2]

        File? InputSpreadSheet
        File? DemuxSampleSheet
        String? PreviousFastqDir
        String? IlluminaDir
        String? XferLabel
        String? DragenEnv

        Boolean DataTransfer
        Boolean RmRunDir

        String OutputDir
        String Queue
        String JobGroup
        String DragenMEM
        String DragenQueue
        String DragenDockerImage

        Int DragenCPU

        String SomaRepo
        String CoverageBed  = SomaRepo + "/accessory_files/SOMA.all.bed"
        String HaplotectBed = SomaRepo + "/accessory_files/SOMA.haplotect.bed"
        String QC_py        = SomaRepo + "/scripts/QC_metrics.py"

        String CovLevels = "100,500,1000,1500"
    }

    String DragenReference = "/storage1/fs1/duncavagee/Active/SEQ/reference/dragen_hg38v4.3.6"
    String Reference       = "/storage1/fs1/duncavagee/Active/SEQ/reference/hg38/sequence/hg38_mgi_patch.fa"
    String ReferenceDict   = "/storage1/fs1/duncavagee/Active/SEQ/reference/hg38/sequence/hg38_mgi_patch.dict"

    String DemuxFastqDir = "/storage1/fs1/gtac-mgi/Active/CLE/assay/SOMA/demux_fastq"

    Int readfamilysize  = 1

    if (defined(DemuxSampleSheet)){
        call dragen_demux {
            input: Dir=IlluminaDir,
                   OutputDir=OutputDir,
                   DemuxFastqDir=DemuxFastqDir,
                   SampleSheet=DemuxSampleSheet,
                   DragenCPU=DragenCPU,
                   DragenMEM=DragenMEM,
                   DragenEnv=DragenEnv,
                   DragenDockerImage=DragenDockerImage,
                   queue=DragenQueue,
                   jobGroup=JobGroup
        }

        call prepare_samples {
            input: SampleSheet=SampleSheet,
                   PreviousFastqDir=PreviousFastqDir,
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
                   DragenCPU=DragenCPU,
                   DragenMEM=DragenMEM,
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

    if (DataTransfer) {
        call data_transfer {
            input: QcAll=batch_qc.QC_all,
            #QcFile= OutputDir + '/' + basename(OutputDir) + '_Genoox.xlsx',
            QcFile=batch_qc.QC_file,
            BatchFastqDir= DemuxFastqDir + '/' + basename(OutputDir),
            InputSpreadSheet=InputSpreadSheet,
            XferLabel=XferLabel,
            queue=Queue,
            jobGroup=JobGroup
        }
    }

    if (defined(DemuxSampleSheet)){
        if (RmRunDir) {
            call remove_rundir {
                input: order_by=gather_files.done,
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
         String? DragenEnv
         File? SampleSheet
         String OutputDir
         String DemuxFastqDir
         String DragenDockerImage
         String DragenMEM
         String jobGroup
         String queue
         Int DragenCPU
     }

     String LocalFastqDir   = "/staging/runs/Soma/demux_fastq/" + basename(OutputDir)
     String OutputFastqDir  = DemuxFastqDir + "/" + basename(OutputDir)
     String OutputReportDir = OutputFastqDir + "/Reports"
     String DemuxReportDir  = OutputDir + "/dragen_demux_reports"

     command <<<
         /bin/mkdir ~{LocalFastqDir} && \
         /opt/dragen/4.3.6/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true --sample-sheet ~{SampleSheet} \
         --bcl-input-directory ~{Dir} --intermediate-results-dir ~{LocalFastqDir} --output-directory ~{OutputFastqDir} && \
         /bin/ls ~{OutputFastqDir}/*_R1_001.fastq.gz > Read1_list.txt && \
         /bin/ls ~{OutputFastqDir}/*_R2_001.fastq.gz > Read2_list.txt && \
         /bin/cp -r ~{OutputReportDir} ~{DemuxReportDir}
     >>>

     runtime {
         docker_image: DragenDockerImage
         cpu: DragenCPU
         memory: DragenMEM
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
         String? PreviousFastqDir
         File SampleSheet
         String Fastq1
         String Fastq2
         String jobGroup
         String queue
     }

     command <<<
         /bin/cp ~{Fastq1} 1.tmp.txt
         /bin/cp ~{Fastq2} 2.tmp.txt
         /usr/bin/perl -e 'use File::Basename; open(R1,"1.tmp.txt"); @r1 = <R1>; \
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
                 my $prev_dir = "~{PreviousFastqDir}";
                 if ($prev_dir) {
                     my ($n) = basename($r1) =~/^(\S+?)_/;
                     my @p_r1 = glob($prev_dir."/$n"."*_R1_001.fastq.gz");
                     my @p_r2 = glob($prev_dir."/$n"."*_R2_001.fastq.gz");
                     unless (@p_r1 and @p_r1 == 1 and @p_r2 and @p_r2 == 1) {
                         die "fail to get previous R1 and or R2 for $n";
                     }
                     my $rc1 = system "cat $p_r1[0] >> $r1";
                     my $rc2 = system "cat $p_r2[0] >> $r2";
                     unless ($rc1 == 0 and $rc2 == 0) {
                         die "R1 and or R2 cat failed for $n";
                     }
                 }
                 print join("\t",@l,$r1,$r2),"\n";
             }
             close SS;' > sample_sheet.txt
     >>>

     runtime {
         docker_image: "docker1(ubuntu:xenial)"
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
         String DragenMEM
         String jobGroup
         String queue

         String? DragenEnv
         Int DragenCPU
         Int readfamilysize
     }

     String batch = basename(OutputDir)
     String LocalAlignDir = "/staging/runs/Soma/align/" + batch + "/" + SubDir
     String DragenOutdir = OutputDir + "/" + SubDir + "/dragen"

     command {
         /bin/mkdir -p ${LocalAlignDir} && \
         /bin/mkdir -p ${DragenOutdir}  && \
         /opt/dragen/4.3.6/bin/dragen -r ${DragenRef} --tumor-fastq1 ${fastq1} --tumor-fastq2 ${fastq2} --RGSM-tumor ${SM} --RGID-tumor ${RG} --RGLB-tumor ${LB} \
         --umi-enable true --umi-library-type=random-simplex --umi-min-supporting-reads ${readfamilysize} --umi-metrics-interval-file ${CoverageBed} \
         --enable-map-align true --enable-sort true --enable-map-align-output true --gc-metrics-enable=true \
         --enable-variant-caller=true --vc-target-bed ${CoverageBed} --vc-enable-umi-solid true --vc-enable-triallelic-filter false \
         --vc-combine-phased-variants-distance 3 --vc-enable-orientation-bias-filter true --vc-skip-germline-tagging true --vc-systematic-noise NONE \
         --qc-coverage-ignore-overlaps=true --qc-coverage-region-1 ${CoverageBed} --qc-coverage-reports-1 cov_report \
         --qc-coverage-region-1-thresholds ${CovLevels} \
         --intermediate-results-dir ${LocalAlignDir} --output-dir ${DragenOutdir} --output-file-prefix ${Name} --output-format CRAM
     }

     runtime {
         docker_image: DragenDockerImage 
         cpu: DragenCPU
         memory: DragenMEM
         dragen_env: DragenEnv
         queue: queue
         job_group: jobGroup
     }

     output {
         File cram = "${DragenOutdir}/${Name}_tumor.cram"
         File crai = "${DragenOutdir}/${Name}_tumor.cram.crai"
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
         docker_image: "docker1(abelhj/haplotect:0.3)"
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
         docker_image: "docker1(ubuntu:xenial)"
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
         String BatchDir
         String QC_py
         String queue
         String jobGroup
         String? InputSpreadSheet
     }
     String batch = basename(BatchDir)

     command {
         if [ -n "$(/bin/ls -d ${BatchDir}/G*)" ]; then
             /bin/chmod -R 666 ${BatchDir}/G*
         fi
         if [ -n "${InputSpreadSheet}" ]; then
             /usr/bin/python3 ${QC_py} -s ${InputSpreadSheet} -d ${BatchDir}
         else
             /usr/bin/python3 ${QC_py} -d ${BatchDir}
         fi
     }
     runtime {
         docker_image: "docker1(catgumag/pandas-scibioxl:20220107)"
         memory: "4 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         File  QC_all  = "${BatchDir}/${batch}_QC.xlsx"
         File? QC_file = "${BatchDir}/${batch}_Genoox.xlsx"
     }
}

task remove_rundir {
     input {
         Array[String] order_by
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
         String? InputSpreadSheet
         String? XferLabel
         String? QcFile
         String QcAll
         String BatchFastqDir
         String queue
         String jobGroup
     }
     
     command {
         set -eo pipefail && \
         /bin/mkdir xfer_staging && \
         if [ -n "${InputSpreadSheet}" ]; then
             /bin/cp ${QcFile} xfer_staging
         fi
         /bin/cp ${BatchFastqDir}/*.fastq.gz xfer_staging && \
         /usr/local/bin/aws s3 cp xfer_staging s3://genoox-upload-wustl/gtacmgi/${XferLabel} --exclude "Undetermined*" --recursive && \
         /usr/bin/touch done.txt && \
         /usr/local/bin/aws s3 cp done.txt s3://genoox-upload-wustl/gtacmgi/${XferLabel}
     }
     runtime {
         docker_image: "docker1(mgibio/data-transfer-aws:v1)"
         memory: "8 G"
         queue: queue
         job_group: jobGroup
     }
     output {
         String done = stdout()
     }
}
