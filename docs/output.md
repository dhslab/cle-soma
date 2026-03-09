# dhslab/cle-soma: Output

All paths below are relative to `--outdir`.

## Per-sample outputs

- `<sample>_<index>/dragen/`
  - `*_tumor.cram`
  - `*_tumor.cram.crai`
  - DRAGEN metrics and logs
- `<sample>_<index>/`
  - `*.haplotect.txt`
  - `*.haplotectloci.txt`

## Batch outputs

- `<batch>_QC.xlsx`
- `<batch>_Genoox.xlsx` (optional; generated when `--input_spreadsheet` is provided)

## Optional transfer artifacts

When `--data_transfer true`:
- Fastq files from `${params.demux_outdir}/<batch>` are staged and uploaded to S3.
- `done.txt` is uploaded as transfer completion marker.

## Notes

- If `--demux_samplesheet` is set, demultiplexed reads are produced before alignment.
- If `--rm_rundir true` and demultiplex mode is enabled, the Illumina run directory is removed at the end.
