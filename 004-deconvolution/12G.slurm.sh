#!/bin/bash
#SBATCH --job-name=12G
#SBATCH --output=12G.out
#SBATCH --error=12G.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

module load cellranger/7.2.0
cellranger count --id=12G --transcriptome=resources/refdata-gex-GRCh38-2020-A/ --fastqs=resources/fastqs/12G --sample=12G --localcores=6 --localmem=25 && rm -rf results/cellranger_count/pool12G && mv pool12G results/cellranger_count/
