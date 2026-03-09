process REMOVE_RUNDIR {
    label 'process_single'
    container params.ubuntu_container

    input:
    val(rundir)

    output:
    path('versions.yml'), emit: versions

    script:
    """
    if [ -d "${rundir}" ]; then
        rm -rf "${rundir}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rm: system
    END_VERSIONS
    """
}
