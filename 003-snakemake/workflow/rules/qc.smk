rule qc_initial:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        histogram="figures/initial_qc/histogram/{sample}.png",
        spatial="figures/initial_qc/spatial/metrics/{sample}.png",
        features_mean_top_100="results/initial_qc/features_mean_top_100/{sample}.tsv",
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
        histogram="figures/filtered_genes_qc/histogram/{sample}.png",
        spatial="figures/filtered_genes_qc/spatial/metrics/{sample}.png",
        features_mean_top_100="results/filtered_genes_qc/features_mean_top_100/{sample}.tsv",
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
        expand("results/initial_qc/features_mean_top_100/{sample}.tsv", sample = samples.index.tolist()),
    output:
        tsv="results/initial_qc/features_mean_top_100/_detection_rate.tsv",
    log:
        "logs/qc/most_abundant_feature_detection_rate.log",
    script:
        "../../scripts/most_abundant_feature_detection_rate.py"

rule most_detected_most_abundant_features:
    input:
        tsv="results/initial_qc/features_mean_top_100/_detection_rate.tsv",
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
