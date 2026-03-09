process DRAGEN_DEMULTIPLEX {
    label 'dragen'
    container params.dragen_container

    input:
    path(demux_samplesheet)
    path(illumina_rundir)

    output:
    path('Read1_list.txt'), emit: read1
    path('Read2_list.txt'), emit: read2
    path('Reports')       , emit: reports
    path('versions.yml')  , emit: versions

    script:
    def batch_name = params.batch_name ?: params.outdir.tokenize('/').last()
    def local_fastq_dir = "/staging/runs/Soma/demux_fastq/${batch_name}"
    def output_fastq_dir = "${params.demux_outdir}/${batch_name}"
    """
    mkdir -p ${local_fastq_dir}
    /opt/dragen/4.3.6/bin/dragen \
        --bcl-conversion-only true \
        --bcl-only-matched-reads true \
        --strict-mode true \
        --sample-sheet ${demux_samplesheet} \
        --bcl-input-directory ${illumina_rundir} \
        --intermediate-results-dir ${local_fastq_dir} \
        --output-directory ${output_fastq_dir}

    ls ${output_fastq_dir}/*_R1_001.fastq.gz > Read1_list.txt
    ls ${output_fastq_dir}/*_R2_001.fastq.gz > Read2_list.txt
    cp -r ${output_fastq_dir}/Reports Reports

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/dragen/4.3.6/bin/dragen --version | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}
