process RUN_HAPLOTECT {
    tag "${meta.id}"
    label 'process_single'
    container params.haplotect_container

    input:
    tuple val(meta), path(cram), path(crai)
    path(haplotect_bed)
    path(reference)
    path(reference_dict)

    output:
    tuple val(meta), path("${meta.id}.haplotect.txt"), path("${meta.id}.haplotectloci.txt"), emit: haplotect_files
    path('versions.yml'), emit: versions

    script:
    """
    awk -v OFS="\\t" '{ \$2=\$2-1; print; }' ${haplotect_bed} > pos.bed

    /usr/local/openjdk-8/bin/java -Xmx6g \
        -jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar Haplotect \
        -I ${cram} \
        -R ${reference} \
        --sequence-dictionary ${reference_dict} \
        -mmq 20 \
        -mbq 20 \
        -max-depth-per-sample 10000 \
        -gstol 0.001 \
        -mr 10 \
        -htp ${haplotect_bed} \
        -L pos.bed \
        -outPrefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
