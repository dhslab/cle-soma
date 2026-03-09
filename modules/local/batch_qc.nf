process BATCH_QC {
    label 'process_single'
    container params.pandas_container

    input:
    val(order_by)
    path(qc_script)
    val(input_spreadsheet)

    output:
    path("*_QC.xlsx"), emit: qc_all
    path("*_Genoox.xlsx"), emit: qc_file, optional: true
    path('versions.yml'), emit: versions

    script:
    def batch_name = params.batch_name ?: params.outdir.tokenize('/').last()
    def input_spreadsheet_opt = input_spreadsheet ? "-s ${input_spreadsheet}" : ""
    """
    if [ -n "\$(ls -d ${params.outdir}/G* 2>/dev/null)" ]; then
        chmod 666 ${params.outdir}/G*/*.txt || true
    fi

    python3 ${qc_script} -d ${params.outdir} ${input_spreadsheet_opt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
