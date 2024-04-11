## Overview

* Run [Cell Ranger][cellranger-github] to prepare the BAM file for Cellsnp-lite.
* Run [Cellsnp-lite][cellsnp-lite-github] to prepare input data for Vireo.
* Run [Vireo][vireo-github] for Demultiplexing pooled scRNA-seq data with or without genotype reference.

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

## Snakemake

### Environment

```bash
miniforge_activate_base # alias
conda activate envs/snakemake
```

### Dry-run

```bash
snakemake --dry-run --profile slurm
```

### Running

```bash
snakemake --use-envmodules --profile slurm
```

[cellranger-github]: https://github.com/10XGenomics/cellranger
[cellsnp-lite-github]: https://github.com/single-cell-genetics/cellsnp-lite
[vireo-github]: https://github.com/single-cell-genetics/vireo
