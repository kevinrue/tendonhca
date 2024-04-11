#!/bin/bash
#SBATCH --job-name=cellsnp-lite
#SBATCH --output=results/cellsnp-lite/12G.out
#SBATCH --error=results/cellsnp-lite/12G.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=1:00:00

# environment
eval "$(/project/sims-lab/albrecht/miniforge3/bin/conda shell.bash hook)" && conda activate base
conda activate cellsnp_lite

# input
BAM=results/cellranger_count/12G/outs/possorted_genome_bam.bam
BARCODE=results/cellranger_count/12G/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
OUT_DIR=results/cellsnp-lite/12G

# cmd
cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --gzip
