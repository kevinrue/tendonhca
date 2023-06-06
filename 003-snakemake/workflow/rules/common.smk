import yaml
import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_final_output():
    final_output = ["results/spaceranger_stats/runtime.tsv", # spaceranger_runtime
                    "results/spaceranger_stats/total_counts.tsv", # spaceranger_total_counts
                    "results/initial_qc/features_mean_top_100/_detection_rate.tsv", # most_abundant_feature_detection_rate
                    "annotations/genes_mitochondrial.tsv", # genes_mitochondrial
                    "annotations/genes_ribosomal.tsv", # genes_ribosomal
                    "figures/spatial/curated_celltype_markers_counts_full",
                    "figures/spatial/curated_celltype_markers_counts_quantile",
                    "figures/spatial/curated_celltype_markers_counts_log1p"]
    # initial quality control
    final_output.append(expand(
        "results/initial_qc/features_mean_top_100/{sample}.tsv", # qc_initial
        sample=samples.index.tolist(),
    ))
    final_output.append(expand(
        "figures/initial_qc/spatial/slide/{sample}.png", # spatial_basic_qc
        sample=samples.index.tolist(),
    ))
    # quality control after filtering genes
    final_output.append(expand(
        "figures/qc_filtered_genes/histogram/{sample}.png",
        sample=samples.index.tolist(),
    ))
    # dimensionality reduction
    final_output.append(expand(
        "figures/filtered_genes/{sample}/dimred/pca.png",
        sample=samples.index.tolist(),
    ))
    final_output.append(expand(
        "figures/filtered_genes/{sample}/dimred/umap.png",
        sample=samples.index.tolist(),
    ))
    final_output.append(expand(
        "figures/filtered_genes/{sample}/dimred/spatial_clusters_image.png",
        sample=samples.index.tolist(),
    ))
    # geneset enrichment
    final_output.append(expand(
        "figures/filtered_genes/{sample}/dimred/spatial_clusters_labelled_image.png",
        sample=samples.index.tolist(),
    ))
    # most abundant features
    final_output.append(expand(
        "results/qc_filtered_genes/features_mean_top_100/{sample}.tsv",
        sample=samples.index.tolist(),
    ))
    final_output.append(expand(
        "figures/spatial/most_detected_most_abundant_features/{sample}_image.png", # spatial_most_detected
        sample=samples.index.tolist(),
    ))
    return final_output

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}

def get_image(wildcards):
    u = samples.loc[wildcards.sample]
    return {"image": u["image"]}

def get_slide(wildcards):
    u = samples.loc[wildcards.sample]
    return u["slide"]

def get_area(wildcards):
    u = samples.loc[wildcards.sample]
    return u["area"]