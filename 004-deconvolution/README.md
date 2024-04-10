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
snakemake --dry-run
```

### Running

```bash
snakemake --slurm --use-envmodules --use-conda
```