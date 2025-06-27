
## History

### 25 Jun 2025

Plan:

- Follow GATK best practices to [Identify short variants (SNPs and Indels) in RNAseq data](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels).
  + Mapping to the reference with 'STAR'.
  + Data cleanup with 'MergeBamAlignment' and 'MarkDuplicates'.
  + SplitNCigarReads to split reads at N bases into supplementary alignments.
  + Base Quality Recalibration using 'BaseRecalibrator', 'Apply Recalibration', 'AnalyzeCovariates'.
  + Variant Calling using 'HaplotypeCaller'.
  + Variant Filtering using 'VariantFiltration'.

Copy Snakemake template repository from <https://github.com/snakemake-workflows/snakemake-workflow-template>.

```bash
wget https://github.com/snakemake-workflows/snakemake-workflow-template/archive/refs/heads/main.zip
unzip main.zip
rm main.zip
mv snakemake-workflow-template-main/* .
mv snakemake-workflow-template-main/.* .
rmdir snakemake-workflow-template-main
```

#### Genome FASTA file

Edit workflow to download the reference genome from Ensembl.

Note that the 'primary assembly' genome is used, as per recommendations in the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

Using VSCode command palette, create a Conda environment for the project.

Install latest version of Snakemake.

```bash
mamba install -c conda-forge snakemake
pip install snakemake-executor-plugin-slurm
```

Reminder of my Snakemake configuration:

```yaml
executor: slurm
jobs: 100
latency-wait: 15
default-resources:
  - slurm_partition=short
  - runtime=10
  - mem=2G
  - tmpdir=system_tmpdir
rerun-incomplete: True
printshellcmds: True
slurm-keep-successful-logs: True
```

Check.

```bash
snakemake --dry-run
```

Run.

```bash
rm nohup.out && nohup snakemake --sdm conda &
tail -f nohup.out
```

Fix error:

> Conda must be version 24.7.1 or later, found version 23.3.1. Please update conda to the latest version. Note that you can also install conda into the snakemake environment without modifying your main conda installation.

```bash
mamba install -c conda-forge conda
```

Run again.

```bash
rm nohup.out && nohup snakemake --sdm conda &
tail -f nohup.out
```

Fix warning:

> Your conda installation is not configured to use strict channel priorities. This is however important for having robust and correct environments (for details, see https://conda-forge.org/docs/user/tipsandtricks.html). Please consider to configure strict priorities by executing 'conda config --set channel_priority strict'.

```bash
conda config --set channel_priority strict
```

#### Genome GTF file

Extend workflow to download the genome annotation from Ensembl.

Note that the 'patch hapl scaff' annotations are used, as per recommendations in the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
Specifically, they recommend using the most comprehensive annotation available, and I thought it was better to have annotations for sequences that are not part of the primary assembly, rather than missing annotations for sequences that are part of the primary assembly.

Run workflow again to download the genome annotation.

#### Index genome with STAR

Extend workflow to index the genome using STAR.

Note that STAR needs to know the read length for the `--sjdbOverhang` option.
I don't have read access to the FASTA files, so I provisionally assumed 150bp while sending an email to get the number (and access to the files).

#### Map reads to the genome with STAR

Extend workflow to map the reads to the genome using STAR.

Run workflow again to map the reads to the genome.

```bash
rm nohup.out && nohup snakemake --sdm conda &
```

### 26 Jun 2025

#### Mark duplicates

Extend workflow to mark duplicates using GATK's `MarkDuplicates`.

Run workflow again to map the reads to the genome.

```bash
rm nohup.out && nohup snakemake --sdm conda &
```

Woops.
I made STAR output BAM files, but called those SAM files in the workflow.
Renaming the output files to BAM files triggered Snakemake warnings forcing me to rerun the whole workflow from the top.
Weird, STAR seems to take longer to run than yesterday (over an hour).
Or maybe it's because it is making BAM files instead of SAM files, which it happened while I was walking home from work.

Woops #2.
The MarkDuplicates rule failed because the BAM files do not contain read group information.
That works out well since I was working on a rule to convert the FASTQ files to uBAM files with read group information.

#### Make uBAM files

Based on How-to <https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM>
and meaning of read groups <https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups>.

I tweaked `samples.tsv` to include the read group information, and added a rule to convert the FASTQ files to uBAM files using GATK's `FastqToSam`.

```bash
rm nohup.out && nohup snakemake --sdm conda &
```

#### Create sequence dictionary to use with GATK tools

Extend workflow to create a sequence dictionary for the genome FASTA file using GATK's `CreateSequenceDictionary`.

I'm still waiting for the latest run of the workflow to finish, so I can run this rule.

#### Split reads at N bases

Extend workflow to split reads at N bases into supplementary alignments using GATK's `SplitNCigarReads`.

I'm still waiting for the latest run of the workflow to finish, so I can re-run the workflow after fixing the creation of a sequence dictionary.

#### Base Quality Recalibration

Extend workflow to perform Base Quality Recalibration using GATK's `BaseRecalibrator` and `ApplyBQSR`.

I'm still waiting for the latest run of the workflow to finish, so I can re-run the workflow after fixing the creation of a sequence dictionary.

### 27 Jun 2025

#### Re-run workflow after fix

Re-running the workflow after fixing the creation of a sequence dictionary.

Re-run includes:

- Marking duplicates

#### TODO

Extended the pipeline up to the merging of filtered variant calls.

## Resources

- GATK best practices for RNAseq: <https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels>
- STAR manual <https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf>
- Minor allele frequency for 1000genomes without downloading GBs of VCF files? <https://bioconductor.org/packages/release/data/annotation/html/MafDb.1Kgenomes.phase3.GRCh38.html>
