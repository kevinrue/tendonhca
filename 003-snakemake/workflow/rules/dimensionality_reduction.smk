rule dimred_filtered_genes_counts:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        pca="figures/filtered_genes/{sample}/dimred/pca.png",
        umap="figures/filtered_genes/{sample}/dimred/umap.png",
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
