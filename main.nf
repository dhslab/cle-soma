#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dhslab/cle-soma
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WDL to Nextflow DSL2 conversion
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { SOMA } from './workflows/soma'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_soma_pipeline/main'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_soma_pipeline/main'

workflow {
    ch_versions = Channel.empty()

    PIPELINE_INITIALISATION(
        params.help,
        params.version,
        params.validate_params,
        params.input,
        params.outdir
    )
    ch_versions = ch_versions.mix(PIPELINE_INITIALISATION.out.versions)

    SOMA(
        PIPELINE_INITIALISATION.out.input
    )
    ch_versions = ch_versions.mix(SOMA.out.versions)

    PIPELINE_COMPLETION(
        ch_versions
    )
}
