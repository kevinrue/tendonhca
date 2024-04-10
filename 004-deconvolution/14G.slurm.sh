#!/bin/bash
#SBATCH --job-name=14G
#SBATCH --output=14G.out
#SBATCH --error=14G.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

module load cellranger/7.2.0
cellranger count --id=14G --transcriptome=resources/refdata-gex-GRCh38-2020-A/ --fastqs=resources/fastqs/14G --sample=14G --localcores=6 --localmem=25 && rm -rf results/cellranger_count/pool14G && mv pool14G results/cellranger_count/
