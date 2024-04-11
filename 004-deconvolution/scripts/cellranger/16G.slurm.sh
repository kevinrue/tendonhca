#!/bin/bash
#SBATCH --job-name=16G
#SBATCH --output=16G.out
#SBATCH --error=16G.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

module load cellranger/7.2.0
cellranger count --id=16G --transcriptome=resources/refdata-gex-GRCh38-2020-A/ --fastqs=resources/fastqs/16G --sample=16G --localcores=6 --localmem=25 && rm -rf results/cellranger_count/16G && mv 16G results/cellranger_count/
