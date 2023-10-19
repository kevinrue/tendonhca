# General activation workflow

```bash
miniforge_activate_base
cd /ceph/project/tendonhca/albrecht/003-snakemake
conda activate envs/cell2location-jaxcuda
sbatch notebooks/hamstring-gpu-slurm.sh
```

# Attempt N

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - pip
  - python
  - anndata
  - leidenalg
  - scanpy
  - squidpy  
  - pip:
    - ipykernel # needed for jupyter notebook
    - nbconvert
    - jax[cuda12_pip]
    - cell2location
```


```bash
conda activate envs/cell2location-jaxcuda
sbatch notebooks/hamstring-gpu-slurm.sh
```

See also:

- <https://cell2location.readthedocs.io/en/latest/commonerrors.html?highlight=cuda#training-cell2location-on-gpu-takes-forever-50-hours>
