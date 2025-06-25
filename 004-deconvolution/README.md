## Overview

* Run [Cell Ranger][cellranger-github] to prepare the BAM file for Cellsnp-lite.
* Run [Cellsnp-lite][cellsnp-lite-github] to prepare input data for Vireo.
  * [Mode 2a: droplet based single cells without given SNPs][cellsnp-lite-mode2a]
* Run [Vireo][vireo-github] for Demultiplexing pooled scRNA-seq data with or without genotype reference.

Other ideas:

- <https://github.com/jon-xu/scSplit>

## Snakemake

### Environment

```bash
miniforge_activate_base # alias
conda activate envs/snakemake
```

### Dry-run

```bash
snakemake --dry-run
```

### Running

```bash
snakemake --use-envmodules
```

Workaround to avoid issue #1 (below):

```bash
module load cellranger/7.2.0
snakemake
```

## Conda

Activate Miniforge.

```bash
miniforge_activate_base
```

Create the Conda environments.

```bash
mamba env create -f conda.yml
mamba env create -f cellsnp_lite.yml
```

Activate the `vireo` environment.

```bash
conda activate vireo
```

Activate the `cellsnp_lite` environment.

```bash
conda activate cellsnp_lite
```

## Notes

### Issue #1: Snakemake and Slurm

See:

* https://github.com/snakemake/snakemake/issues/2802
* https://github.com/snakemake/snakemake-executor-plugin-slurm/issues/63

[cellranger-github]: https://github.com/10XGenomics/cellranger
[cellsnp-lite-github]: https://github.com/single-cell-genetics/cellsnp-lite
[vireo-github]: https://github.com/single-cell-genetics/vireo
[cellsnp-lite-mode2a]: https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html#mode-2a-droplet-based-single-cells-without-given-snps