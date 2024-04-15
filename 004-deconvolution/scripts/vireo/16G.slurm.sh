#!/bin/bash
#SBATCH --job-name=vireo
#SBATCH --output=results/vireo/16G.out
#SBATCH --error=results/vireo/16G.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=23:00:00

# environment
eval "$(/project/sims-lab/albrecht/miniforge3/bin/conda shell.bash hook)" && conda activate base
conda activate vireo

# input
CELL_DATA=results/cellsnp-lite/16G
n_donor=7
OUT_DIR=results/vireo/16G

# cmd
vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR
