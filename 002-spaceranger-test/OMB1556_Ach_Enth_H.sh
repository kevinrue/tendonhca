#!/bin/sh

#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=kevin.rue-albrecht@imm.ox.ac.uk
#SBATCH --mail-type=begin,end,fail
#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=6
#SBATCH --mem=25G

cd /project/tendonhca/albrecht/002-spaceranger-test

module load spaceranger/1.3.1

spaceranger count --id=OMB1556_Ach_Enth_H \
  --transcriptome=/project/tendonhca/shared/spatial/analysis/refdata-gex-GRCh38-2020-A/ \
  --fastqs=./data.dir/OMB1556_Ach_Enth/ \
  --sample=OMB1556_Ach_Enth_H \
  --image=/project/tendonhca/shared/spatial/analysis/achilles/images/OMB1556_Ach_Enth.jpeg \
  --slide=V12J03-133 \
  --area=C1 \
  --localcores=6 \
  --localmem=25
