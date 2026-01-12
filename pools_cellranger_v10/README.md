## Note

Repeat of `006-single-nuc-cellranger/` using [snakemake-workflows/cellranger-multi](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-multi) branch `fix/account-for-different-output-file-paths-with-cellranger-10-0-0` and CellRanger v10.

## History

### 12 Jan 2026

```bash
snakedeploy deploy-workflow https://github.com/snakemake-workflows/cellranger-multi . --branch fix/account-for-different-output-file-paths-with-cellranger-10-0-0
```

Updated config files:

`config/pools.tsv`

- samples and file paths

`config/config.yaml`

- path to `refdata-gex-GRCh38-2020-A/`

```bash
export CELLRANGER_TARBALL="/ceph/project/tendonhca/albrecht/resources/cellranger-10.0.0.tar.gz"
snakemake --dry-run
nohup snakemake --cores all --sdm conda apptainer --default-resources runtime=720 &
```
