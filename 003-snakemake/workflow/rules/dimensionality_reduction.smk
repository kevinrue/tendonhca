rule dimred_filtered_genes_counts:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        genes='annotations/genes.tsv',
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        pca="figures/filtered_genes/{sample}/dimred/pca.png",
        umap="figures/filtered_genes/{sample}/dimred/umap.png",
        spatial_clusters="figures/filtered_genes/{sample}/dimred/spatial_clusters.png",
        highly_variable_genes="figures/filtered_genes/{sample}/dimred/highly_variable_genes.png",
        pca_variance_ratio="figures/filtered_genes/{sample}/dimred/pca_variance_ratio.png",
        rank_genes_groups="results/filtered_genes/{sample}/markers/rank_genes_groups.tsv.gz",
    params:
        samples=config["samples"],
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        "logs/dimred_filtered_genes_counts/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    script:
        "../../scripts/dimred_filtered_genes.py"
