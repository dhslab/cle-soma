process DRAGEN_ALIGNMENT {
    tag "${meta.id}"
    label 'dragen'
    container params.dragen_container

    input:
    tuple val(meta), path(fastq1), path(fastq2)
    path(coverage_bed)

    output:
    tuple val(meta), path("dragen/${meta.id}_tumor.cram"), path("dragen/${meta.id}_tumor.cram.crai"), emit: cram_files
    tuple val(meta), path("dragen/*"), emit: dragen_files
    path('versions.yml'), emit: versions

    script:
    def rg = "${meta.flowcell}.${meta.lane}.${meta.index}"
    def batch = params.batch_name ?: params.outdir.tokenize('/').last()
    def local_align_dir = "/staging/runs/Soma/align/${batch}/${meta.subdir}"

    """
    mkdir -p ${local_align_dir}
    mkdir -p dragen

    /opt/dragen/4.3.6/bin/dragen \
        -r ${params.dragen_reference} \
        --tumor-fastq1 ${fastq1} \
        --tumor-fastq2 ${fastq2} \
        --RGSM-tumor ${meta.sm} \
        --RGID-tumor ${rg} \
        --RGLB-tumor ${meta.lb}.${meta.index} \
        --umi-enable true \
        --umi-library-type random-simplex \
        --umi-min-supporting-reads ${params.read_family_size} \
        --umi-metrics-interval-file ${coverage_bed} \
        --enable-map-align true \
        --enable-sort true \
        --enable-map-align-output true \
        --gc-metrics-enable=true \
        --enable-variant-caller=true \
        --vc-target-bed ${coverage_bed} \
        --vc-enable-umi-solid true \
        --vc-enable-triallelic-filter false \
        --vc-combine-phased-variants-distance 3 \
        --vc-enable-orientation-bias-filter true \
        --vc-skip-germline-tagging true \
        --vc-systematic-noise NONE \
        --qc-coverage-ignore-overlaps=true \
        --qc-coverage-region-1 ${coverage_bed} \
        --qc-coverage-reports-1 cov_report \
        --qc-coverage-region-1-thresholds ${params.cov_levels} \
        --intermediate-results-dir ${local_align_dir} \
        --output-dir dragen \
        --output-file-prefix ${meta.id} \
        --output-format CRAM

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/dragen/4.3.6/bin/dragen --version | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}
