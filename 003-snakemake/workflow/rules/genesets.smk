rule markers_fgsea:
    input:
        rank_genes_groups="results/filtered_genes/{sample}/markers/rank_genes_groups.tsv.gz",
    output:
        fgsea="results/filtered_genes/{sample}/markers/fgsea.tsv",
    log:
        "logs/genesets/markers_fgsea_{sample}.log",
    threads: 6
    resources:
        mem_mb=4 * 1024,
        runtime="30m",
    conda:
        "../../envs/r-env.yaml",
    script:
        "../../scripts/markers_fgsea.R"
