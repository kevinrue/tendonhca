rule basic_qc:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        total_counts_n_genes_by_counts="results/qc/total_counts_n_genes_by_counts/{sample}.png",
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
