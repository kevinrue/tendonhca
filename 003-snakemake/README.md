
## Usage

```bash
snakemake --slurm --use-envmodules --use-conda
```

Namely:

- The workload manager available on the HPC cluster used to run this workflow is SLURM.
  - DRMAA is available on the HPC cluster used to run this workflow,
    but fails to pass the requested walltime for `spaceranger count`,
    resulting in jobs being killed.
- Spacer Ranger is available as a module on the HPC cluster used to run this workflow.
- Conda (more specifically Mamba) environments are used for the following dependencies:
  - Scanpy
