rule markers_fgsea:
    input:
        rank_genes_groups="results/filtered_genes/{sample}/markers/rank_genes_groups.tsv.gz",
    output:
        fgsea="results/filtered_genes/{sample}/markers/fgsea.tsv",
    log:
        "logs/genesets/markers_fgsea_{sample}.log",
    threads: 1
    resources:
        mem_mb=4 * 1024,
        runtime=30,
    conda:
        "../../envs/r-env.yaml",
    script:
        "../../scripts/markers_fgsea.R"

rule spatial_clusters_labelled:
    input:
        fgsea="results/filtered_genes/{sample}/markers/fgsea.tsv",
        genes='annotations/genes.tsv',
        mitochondrial='annotations/genes_mitochondrial.tsv',
        ribosomal='annotations/genes_ribosomal.tsv',
    output:
        png="figures/filtered_genes/{sample}/dimred/spatial_clusters_labelled.png",
    log:
        "logs/spatial_clusters_labelled/{sample}.log",
    threads: 1
    resources:
        mem_mb=4 * 1024,
        runtime=10,
    conda:
        "../../envs/scanpy-env.yaml",
    script:
        "../../scripts/autolabeller.py"
