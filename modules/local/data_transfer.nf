process DATA_TRANSFER {
    label 'process_single'
    container params.transfer_container

    input:
    path(qc_all)
    path(qc_file)
    path(fastqs)
    val(input_spreadsheet)

    output:
    path('done.txt'), emit: done
    path('versions.yml'), emit: versions

    script:
    def batch_name = params.batch_name ?: params.outdir.tokenize('/').last()
    """
    set -euo pipefail
    mkdir xfer_staging

    # qc_file is optional (QC Metrics - qPCR / Genoox.xlsx)
    if [ -n "${qc_file}" ]; then
        cp ${qc_file} xfer_staging/
    fi

    cp ${fastqs} xfer_staging/
    aws s3 cp xfer_staging s3://genoox-upload-wustl/gtacmgi/${params.xfer_label} --exclude "Undetermined*" --recursive
    touch done.txt
    aws s3 cp done.txt s3://genoox-upload-wustl/gtacmgi/${params.xfer_label}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(aws --version 2>&1)
    END_VERSIONS
    """
}
