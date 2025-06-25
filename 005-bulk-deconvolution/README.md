
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

Edit workflow to download the reference genome  from Ensembl.

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
