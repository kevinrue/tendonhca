## Note

Repeat of `006-single-nuc-cellranger/` using [snakemake-workflows/cellranger-multi](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-multi).

## History

### 09 Jan 2026

```bash
cd /ceph/project/tendonhca/albrecht/resources
curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1768007656&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=IEn7SIu3Tt9KwaYFnBDS1IFQSnJ8Al21EDDwYXSmRI79SRjOUy5bAcr6UvQ4VYKCDQS0O2X0gNztPoyFfGsA3Y2MVj~G1xEmQMFLC~fpV1bcO6B1D3w3y3HbY3r3LlSRW2uXnFLoaH-PJP5zAsqJMuKzZRKQZ3YKrK1DUaZoeiwUTVG8q8b0ab6EgSOIC6vEcT6bvcWvIQ02CVQfXDGSy5DjH0hqdq0CQnq-ISix-l4Ua9ILJMicxVFyxfwbipN-QU1IFc68GuZ-NEx7dpuldvK9yyIc4DgbG1-TFfGa~-aITzBt1VqkeNn1fBeGFw4TTXY2FUPbuV4a4yWa3iUjaQ__"
```

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
snakemake --cores all --sdm conda apptainer
```
