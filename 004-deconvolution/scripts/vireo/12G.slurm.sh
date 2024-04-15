#!/bin/bash
#SBATCH --job-name=vireo
#SBATCH --output=results/vireo/12G.out
#SBATCH --error=results/vireo/12G.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

# environment
eval "$(/project/sims-lab/albrecht/miniforge3/bin/conda shell.bash hook)" && conda activate base
conda activate vireo

# input
CELL_DATA=results/cellsnp-lite/12G
n_donor=8
OUT_DIR=results/vireo/12G

# cmd
vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR
