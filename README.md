# cle-soma

SOMA assay pipeline in Nextflow DSL2.

## Repository layout

- `main.nf`, `nextflow.config`: Nextflow entrypoint and configuration
- `workflows/`, `subworkflows/`, `modules/`: DSL2 workflow implementation
- `assets/data/`: assay BED files and static pipeline data
- `docs/`: Usage and output documentation
- `legacy/wdl/`: Archived original WDL implementation artifacts

## Start

```bash
nextflow run . -profile docker,ris --input /path/to/samples.tsv --outdir /path/to/output
```

See:
- `docs/usage.md`
- `docs/output.md`
