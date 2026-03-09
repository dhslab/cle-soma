# dhslab/cle-soma: Usage

## Quick Start

```bash
nextflow run /Users/fdu/git/cle-soma \
  -profile docker,ris \
  --input /path/to/soma_samples.tsv \
  --outdir /path/to/output_batch
```

## Required Parameters

- `--input`: Tab-delimited sample sheet with columns:
  - `index  name  RG_ID  RG_FLOWCELL  RG_LANE  RG_LIB  RG_SAMPLE  R1  R2`
- `--outdir`: Batch output directory.

## Optional Parameters

- `--demux_samplesheet`: Enable demultiplex-first workflow.
- `--illumina_rundir`: Illumina run directory; required when `--demux_samplesheet` is used.
- `--demux_outdir`: Where demultiplexed fastqs are written.
- `--input_spreadsheet`: Spreadsheet used by QC script and Genoox export.
- `--data_transfer`: Upload fastq/QC artifacts to S3.
- `--rm_rundir`: Remove Illumina run directory at end of run.
- `--xfer_label`: S3 transfer label suffix.

## Compute Parameters

- `--queue`, `--user_group`, `--job_group_name`
- DRAGEN resources are selected by profile (`dragen2`, `dragen3`, `dragen4`, `dragen5`).

## Reference Parameters

- `--dragen_reference`
- `--reference`
- `--reference_dict`
- `--coverage_bed`
- `--haplotect_bed`
- `--qc_script`
- `--cov_levels`
- `--read_family_size`

## Profiles

- `standard`: local executor
- `docker`: Docker runtime
- `singularity`: Singularity runtime
- `apptainer`: Apptainer runtime
- `ris`: LSF execution profile
- `dragen2`, `dragen3`, `dragen4`, `dragen5`: DRAGEN queue resource presets

Common production invocation:

```bash
nextflow run /Users/fdu/git/cle-soma \
  -profile docker,ris,dragen3 \
  --input /path/to/soma_samples.tsv \
  --outdir /path/to/output_batch
```
