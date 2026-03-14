/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEMULTIPLEX    } from '../subworkflows/local/demultiplex'
include { DRAGEN_ALIGNMENT } from '../modules/local/dragen_alignment'
include { RUN_HAPLOTECT   } from '../modules/local/run_haplotect'
include { GATHER_SAMPLE_FILES } from '../modules/local/gather_sample_files'
include { BATCH_QC        } from '../modules/local/batch_qc'
include { DATA_TRANSFER   } from '../modules/local/data_transfer'
include { REMOVE_RUNDIR   } from '../modules/local/remove_rundir'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SOMA {

    take:
    ch_input_samplesheet

    main:
    ch_versions = Channel.empty()

    if (params.demux_samplesheet && !params.illumina_rundir) {
        error("Please provide --illumina_rundir when using --demux_samplesheet.")
    }

    ch_demux_samplesheet = params.demux_samplesheet
        ? Channel.fromPath(params.demux_samplesheet, checkIfExists: true)
        : Channel.empty()
    ch_illumina_rundir = params.illumina_rundir
        ? Channel.fromPath(params.illumina_rundir, type: 'dir', checkIfExists: true)
        : Channel.empty()
    ch_coverage_bed = Channel.fromPath(params.coverage_bed, checkIfExists: true).collect()
    ch_haplotect_bed = Channel.fromPath(params.haplotect_bed, checkIfExists: true).collect()
    ch_reference = Channel.fromPath(params.reference, checkIfExists: true).collect()
    ch_reference_dict = Channel.fromPath(params.reference_dict, checkIfExists: true).collect()
    ch_qc_script = Channel.fromPath(params.qc_script, checkIfExists: true).collect()

    def chosen_samplesheet = ch_input_samplesheet

    if (params.demux_samplesheet) {
        DEMULTIPLEX(
            ch_input_samplesheet,
            ch_demux_samplesheet,
            ch_illumina_rundir
        )
        ch_versions = ch_versions.mix(DEMULTIPLEX.out.versions)
        chosen_samplesheet = DEMULTIPLEX.out.prepared_samplesheet
    }

    ch_samples = chosen_samplesheet
        .splitCsv(header: false, sep: '\t')
        .map { row ->
            if (row.size() < 9) {
                error("Input samplesheet requires at least 9 tab-delimited columns: index name RG_ID RG_FLOWCELL RG_LANE RG_LIB RG_SAMPLE R1 R2")
            }
            def meta = [
                id     : row[1],
                index  : row[0],
                flowcell: row[3],
                lane   : row[4],
                lb     : row[5],
                sm     : row[6],
                subdir : "${row[1]}_${row[0]}"
            ]
            [meta, file(row[7], checkIfExists: true), file(row[8], checkIfExists: true)]
        }

    DRAGEN_ALIGNMENT(
        ch_samples,
        ch_coverage_bed
    )
    ch_versions = ch_versions.mix(DRAGEN_ALIGNMENT.out.versions)

    RUN_HAPLOTECT(
        DRAGEN_ALIGNMENT.out.cram_files,
        ch_haplotect_bed,
        ch_reference,
        ch_reference_dict
    )
    ch_versions = ch_versions.mix(RUN_HAPLOTECT.out.versions)

    GATHER_SAMPLE_FILES(
        DRAGEN_ALIGNMENT.out.dragen_files
            .map { meta, dragen_files -> [meta.id, meta, dragen_files] }
            .join(
                RUN_HAPLOTECT.out.haplotect_files.map { meta, haplotect_file, haplotect_loci -> [meta.id, meta, haplotect_file, haplotect_loci] }
            )
            .map { sample_id, meta1, dragen_files, meta2, haplotect_file, haplotect_loci -> [meta1, dragen_files, haplotect_file, haplotect_loci] }
    )
    ch_versions = ch_versions.mix(GATHER_SAMPLE_FILES.out.versions)

    BATCH_QC(
        GATHER_SAMPLE_FILES.out.done.collect(),
        ch_qc_script,
        params.input_spreadsheet ?: ""
    )
    ch_versions = ch_versions.mix(BATCH_QC.out.versions)

    if (params.data_transfer) {
        ch_fastqs = ch_samples.map { meta, r1, r2 -> [r1, r2] }.flatten().collect()
        DATA_TRANSFER(
            BATCH_QC.out.qc_all,
            BATCH_QC.out.qc_file.collect(),
            ch_fastqs,
            params.input_spreadsheet ?: ""
        )
        ch_versions = ch_versions.mix(DATA_TRANSFER.out.versions)
    }

    if (params.demux_samplesheet && params.rm_rundir && params.illumina_rundir) {
        REMOVE_RUNDIR(params.illumina_rundir, BATCH_QC.out.versions.collect())
        ch_versions = ch_versions.mix(REMOVE_RUNDIR.out.versions)
    }

    emit:
    versions = ch_versions
}
