rule basic_qc:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        histogram="figures/basic_qc/histogram/{sample}.png",
        spatial="figures/basic_qc/spatial/metrics/{sample}.png",
        features_mean_top_100="results/qc/features_mean_top_100/{sample}.tsv",
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/basic_qc/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/basic_qc.py"

rule most_abundant_feature_detection_rate:
    input:
        expand("results/qc/features_mean_top_100/{sample}.tsv", sample = samples.index.tolist()),
    output:
        tsv="results/qc/features_mean_top_100/_detection_rate.tsv",
    log:
        "logs/qc/most_abundant_feature_detection_rate.log",
    script:
        "../../scripts/most_abundant_feature_detection_rate.py"

rule most_detected_most_abundant_features:
    input:
        tsv="results/qc/features_mean_top_100/_detection_rate.tsv",
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
