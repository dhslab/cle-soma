process GATHER_SAMPLE_FILES {
    tag "${meta.id}"
    label 'process_single'
    container params.ubuntu_container

    publishDir "${params.outdir}/${meta.subdir}", mode: params.publish_dir_mode, saveAs: { filename -> filename.contains('versions.yml') ? null : filename }

    input:
    tuple val(meta), path(dragen_files), path(haplotect_file), path(haplotect_loci)

    output:
    tuple val(meta), val('done'), emit: done
    path(haplotect_file)
    path(haplotect_loci)
    path("dragen/*")
    path('versions.yml'), emit: versions

    script:
    """
    mkdir -p dragen
    cp -f ${dragen_files} dragen/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cp: system
    END_VERSIONS
    """
}
