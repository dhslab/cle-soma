/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DRAGEN_DEMULTIPLEX } from '../../modules/local/dragen_demultiplex'
include { PREPARE_SAMPLES    } from '../../modules/local/prepare_samples'

workflow DEMULTIPLEX {

    take:
    ch_samplesheet
    ch_demux_samplesheet
    ch_illumina_rundir

    main:
    ch_versions = Channel.empty()

    DRAGEN_DEMULTIPLEX(
        ch_demux_samplesheet,
        ch_illumina_rundir
    )
    ch_versions = ch_versions.mix(DRAGEN_DEMULTIPLEX.out.versions)

    PREPARE_SAMPLES(
        ch_samplesheet,
        DRAGEN_DEMULTIPLEX.out.read1,
        DRAGEN_DEMULTIPLEX.out.read2
    )
    ch_versions = ch_versions.mix(PREPARE_SAMPLES.out.versions)

    emit:
    prepared_samplesheet = PREPARE_SAMPLES.out.sample_sheet
    versions             = ch_versions
}
