rule basic_qc:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        total_counts_n_genes_by_counts="results/qc/total_counts_n_genes_by_counts/{sample}.png",
        total_counts_n_genes_by_counts_spatial="results/qc/total_counts_n_genes_by_counts_spatial/{sample}.png",
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
    script:
        "../../scripts/most_abundant_feature_detection_rate.py"
