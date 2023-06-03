rule pca_filtered_genes_counts:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        png="figures/pca/filtered_genes/counts/{sample}.png",
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/pca/filtered_genes/counts/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/pca_filtered_genes_counts.py"

rule pca_filtered_genes_counts_log1p:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        png="figures/pca/filtered_genes/counts_log1p/{sample}.png",
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/pca/filtered_genes/counts_log1p/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/pca_filtered_genes_counts_log1p.py"

rule pca_filtered_genes_counts_log1p_scale:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        png="figures/pca/filtered_genes/counts_log1p_scale/{sample}.png",
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/pca/filtered_genes/counts_log1p_scale/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/pca_filtered_genes_counts_log1p_scale.py"
