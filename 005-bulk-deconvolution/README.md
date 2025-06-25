
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

Extend workflow to download the genome annotation from Ensembl.

Note that the 'patch hapl scaff' annotations are used, as per recommendations in the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
Specifically, they recommend using the most comprehensive annotation available, and I thought it was better to have annotations for sequences that are not part of the primary assembly, rather than missing annotations for sequences that are part of the primary assembly.

Run workflow again to download the genome annotation.

Extend workflow to index the genome using STAR.

Note that STAR needs to know the read length for the `--sjdbOverhang` option.
I don't have read access to the FASTA files, so I provisionally assumed 150bp while sending an email to get the number (and access to the files).

## Resources

- GATK best practices for RNAseq: <https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels>
- STAR manual <https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf>
