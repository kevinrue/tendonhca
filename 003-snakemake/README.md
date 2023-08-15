## Usage

```bash
snakemake --slurm --use-envmodules --use-conda
```

Namely:

- The workload manager available on the HPC cluster used to run this workflow is SLURM.
  - DRMAA is available on the HPC cluster used to run this workflow,
    but fails to pass the requested walltime for `spaceranger count`,
    resulting in jobs being killed.
- Spacer Ranger is available as a module on the HPC cluster used to run this workflow.
- Conda (more specifically Mamba) environments are used for the following dependencies:
  - Scanpy

## Workflow

### Core

- `spaceranger count` is run on each individual sample using information supplied in `config.yaml` and `samples.tsv`.

### Annotations

- Gene information extracted from the input GTF file is stored in `annotations/genes.tsv`.
- Mitochondrial gene identifiers are stored in `annotations/genes_mitochondrial.tsv`.

### Quality control

#### Tables

- The runtime of `spaceranger count` for each sample is collated into the file `results/spaceranger_stats/runtime.tsv`.
- The total count of `spaceranger count` for each sample is collated into the file `results/spaceranger_stats/total_counts.tsv`.
- (deprecate) The top 100 most abundant genes for each sample are stored in `results/qc/features_mean_top_100/{sample}.tsv`.
- (deprecate) The fraction of samples in which genes are detected in the top 100 most abundant is stored in `results/qc/features_mean_top_100/_detection_rate.tsv`.

#### Plots

- Histograms of quality control metrics for individual samples are stored in `figures/qc/total_counts_n_genes_by_counts/{sample}.png`
- Spatial views of quality control metrics for individual samples are stored in `figures/qc/total_counts_n_genes_by_counts_spatial/{sample}.png`.
  + The same plots are combined with a view of the stained tissue in `{sample}_slide.png`.
- Spatial views of the most frequently detected most abundant genes are stored in `figures/spatial/most_detected_most_abundant_features/{sample}.png`.
  + The same plots are combined with a view of the stained tissue in `{sample}_slide.png`.

### Manually curated information

- Spatial views of selected marker genes are stored in:
  - `figures/spatial/curated_celltype_markers_counts_full` (colour scale fitted to full range of counts per spot)
  - `figures/spatial/curated_celltype_markers_counts_quantile` (colour scale capped to 95% maximal value of counts per spot)
  - `figures/spatial/curated_celltype_markers_counts_log1p` (colour scale fitted to log1p-transformed counts per spot)
