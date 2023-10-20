#!/bin/bash

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-01:00:00
#SBATCH --ntasks=2
#SBATCH --mem=12G
#SBATCH --partition=short

# activate miniforge
eval "$(/project/sims-lab/albrecht/miniforge3/bin/conda shell.bash hook)" && conda activate base

conda activate envs/snakemake

snakemake --slurm --use-envmodules --use-conda
