process GATHER_SAMPLE_FILES {
    tag "${meta.id}"
    label 'process_single'
    container params.ubuntu_container

    input:
    tuple val(meta), path(dragen_files), path(haplotect_file), path(haplotect_loci)

    output:
    tuple val(meta), val('done'), emit: done
    path('versions.yml'), emit: versions

    script:
    """
    mkdir -p ${params.outdir}/${meta.subdir}/dragen
    cp -f ${haplotect_file} ${params.outdir}/${meta.subdir}/
    cp -f ${haplotect_loci} ${params.outdir}/${meta.subdir}/
    cp -f dragen_files/* ${params.outdir}/${meta.subdir}/dragen/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cp: system
    END_VERSIONS
    """
}
