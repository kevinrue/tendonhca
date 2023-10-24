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

### Hamstring single-nuclei data as the reference dataset

RDS file located at

```
OneDrive/CZI - Tendon Seed Network/Manuscripts/Hamstring paper/R data/RDS files/HAMSTRING_singlets_ambRNA0.2_res0.15.RDS
```

WARNING when loading the object:

```
The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
which was just loaded, will retire in October 2023.
Please refer to R-spatial evolution reports for details, especially
https://r-spatial.org/r/2023/05/15/evolution4.html.
It may be desirable to make the sf package available;
package maintainers should consider adding sf to Suggests:.
The sp package is now running under evolution status 2
     (status 2 uses the sf package in place of rgdal)
```

## Samples

### Names

- OMB0793_Quad_Enth_T
- OMB0793_Quad_Enth2_T
- OMB1250_Quad_MB_H
- OMB1277_SSP_Enth_H
- OMB1530_LHBT_MB_H
- OMB1530_LHBT_MB2_H
- OMB1541_GluMed_Enth_H
- OMB1541_GluMed_MTJ_H
- OMB1556_Ach_Enth_H
- OMB1556_Ach_MB_H
- OMB1556_Ach_MB2_H
- OMB1556_Ach_MTJ_H
- OMB1574_SSP_Enth_T

### Locations

Tendon types:

- Quadriceps tendon:
  Thick, strong tendon that connects the quadriceps muscles to the patella (kneecap)
  Essential for the extension of the knee joint.
- Supraspinatus tendon:
  One of the four muscles that make up the rotator cuff in the shoulder.
  Located in the upper part of the shoulder, above the spine of the scapula (shoulder blade)
  Plays a key role in the movement and stability of the shoulder joint.
- Long head of the biceps brachii tendon:
  Part of the biceps brachii muscle located in the upper arm.
  The biceps brachii is a two-headed muscle, and the long head is one of those heads.
  It has a distinct tendon that plays a significant role in shoulder and elbow movements.
- Gluteus medius tendon:
  One of the three major muscles in the buttocks, commonly referred to as the glutes.
  Situated on the lateral side of the hip and plays a crucial role in hip stability and movement.
  The gluteus medius tendon plays a crucial role in attaching the muscle to the greater trochanter, a bony prominence on the femur (thigh bone).
- Achilles tendon:
  Thick and strong tendon that connects the calf muscles to the heel bone (calcaneus).

Tendon location:

- Enthesis:
  Point where a tendon or ligament attaches to the bone.
- MB
