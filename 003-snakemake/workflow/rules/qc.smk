rule qc_raw:
    input:
        "results/spaceranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        histogram="figures/qc_raw/histogram/{sample}.png",
        spatial="figures/qc_raw/spatial/metrics/{sample}.png",
        features_mean_top_100="results/qc_raw/features_mean_top_100/{sample}.tsv",
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/qc/qc_raw/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="15m",
    script:
        "../../scripts/qc_raw.py"

rule qc_initial:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        histogram="figures/qc_initial/histogram/{sample}.png",
        spatial="figures/qc_initial/spatial/metrics/{sample}.png",
        features_mean_top_100="results/qc_initial/features_mean_top_100/{sample}.tsv",
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/qc/qc_initial/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/qc_initial.py"

rule qc_filtered_genes:
    input:
        spaceranger="results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        histogram="figures/qc_filtered_genes/histogram/{sample}.png",
        spatial="figures/qc_filtered_genes/spatial/metrics/{sample}.png",
        features_mean_top_100="results/qc_filtered_genes/features_mean_top_100/{sample}.tsv",
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/qc/qc_filtered_genes/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/qc_filtered_genes.py"

rule most_abundant_feature_detection_rate:
    input:
        expand("results/qc_initial/features_mean_top_100/{sample}.tsv", sample = samples.index.tolist()),
    output:
        tsv="results/qc_initial/features_mean_top_100/_detection_rate.tsv",
    log:
        "logs/qc/most_abundant_feature_detection_rate.log",
    script:
        "../../scripts/most_abundant_feature_detection_rate.py"

rule most_detected_most_abundant_features:
    input:
        tsv="results/qc_initial/features_mean_top_100/_detection_rate.tsv",
    output:
        expand("figures/spatial/most_detected_most_abundant_features/{sample}.png", sample = samples.index.tolist()),
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/qc/most_detected_most_abundant_feature.log",
    script:
        "../../scripts/most_abundant_feature_spatial_plot.py"
