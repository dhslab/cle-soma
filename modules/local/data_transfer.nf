process DATA_TRANSFER {
    label 'process_single'
    container params.transfer_container

    input:
    path(qc_all)
    val(input_spreadsheet)

    output:
    path('done.txt'), emit: done
    path('versions.yml'), emit: versions

    script:
    def batch_name = params.batch_name ?: params.outdir.tokenize('/').last()
    def batch_fastq_dir = "${params.demux_outdir}/${batch_name}"
    """
    set -euo pipefail
    mkdir xfer_staging

    if [ -n "${input_spreadsheet}" ] && [ -f "${params.outdir}/${batch_name}_Genoox.xlsx" ]; then
        cp "${params.outdir}/${batch_name}_Genoox.xlsx" xfer_staging/
    fi

    cp ${batch_fastq_dir}/*.fastq.gz xfer_staging/
    aws s3 cp xfer_staging s3://genoox-upload-wustl/gtacmgi/${params.xfer_label} --exclude "Undetermined*" --recursive
    touch done.txt
    aws s3 cp done.txt s3://genoox-upload-wustl/gtacmgi/${params.xfer_label}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(aws --version 2>&1)
    END_VERSIONS
    """
}
