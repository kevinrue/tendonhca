## Note

Repeat of `006-single-nuc-cellranger/` using [snakemake-workflows/cellranger-multi](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-multi).

## History

### 09 Jan 2026

```bash
snakedeploy deploy-workflow https://github.com/snakemake-workflows/cellranger-multi . --tag v2.0.0
```

Updated config files:

`config/pools.tsv`

- samples and file paths

`config/config.yaml`

- path to `refdata-gex-GRCh38-2020-A/`

```bash
export CELLRANGER_TARBALL="/ceph/project/tendonhca/albrecht/resources/cellranger-9.0.1.tar.gz"
snakemake --dry-run
nohup snakemake --cores all --sdm conda apptainer --default-resources runtime=720 &
```

## 12 Jan 2026

```bash
for pool in 12G 14G 16G
do
    cp \
        pools_cellranger/results/cellranger/${pool}/outs/per_sample_outs/${pool}/web_summary.html \
        /ceph/project/tendonhca/datashare/albrecht/pools/cellranger_web_summary/${pool}.html
done
```
