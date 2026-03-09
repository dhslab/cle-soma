/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility subworkflow for cle-soma pipeline lifecycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsHelp         } from 'plugin/nf-validation'
include { paramsSummaryLog   } from 'plugin/nf-validation'
include { validateParameters } from 'plugin/nf-validation'

workflow PIPELINE_INITIALISATION {

    take:
    show_help
    show_version
    validate_params
    input
    outdir

    main:
    ch_versions = Channel.empty()

    if (show_version) {
        log.info "${workflow.manifest.name ?: 'dhslab/cle-soma'} ${workflow.manifest.version ?: 'dev'}"
        System.exit(0)
    }

    if (show_help) {
        def workflow_command = "nextflow run ${workflow.projectDir} --input <samplesheet.tsv> --outdir <output_dir> -profile docker,ris"
        log.info paramsHelp(workflow_command, parameters_schema: 'nextflow_schema.json')
        System.exit(0)
    }

    if (!input) {
        error("Missing required parameter: --input")
    }
    if (!outdir) {
        error("Missing required parameter: --outdir")
    }

    if (params.demux_samplesheet && !params.illumina_rundir) {
        error("--illumina_rundir is required when --demux_samplesheet is provided")
    }

    log.info paramsSummaryLog(workflow, parameters_schema: 'nextflow_schema.json')

    if (validate_params) {
        validateParameters(parameters_schema: 'nextflow_schema.json')
    }

    ch_input = Channel.fromPath(input, checkIfExists: true)

    emit:
    input    = ch_input
    versions = ch_versions
}

workflow PIPELINE_COMPLETION {

    take:
    versions

    main:
    versions.collect().set { _all_versions }
    log.info "Pipeline completion subworkflow initialized."
}
