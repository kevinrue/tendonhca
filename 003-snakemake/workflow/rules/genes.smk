rule gene_table:
    input:
        gtf=config['transcriptome'] + 'genes/genes.gtf'
    output:
        tsv='annotations/genes.tsv'
    script:
        "../../scripts/gtf2tsv.py"

rule genes_mitochondrial:
    input:
        tsv='annotations/genes.tsv'
    output:
        tsv='annotations/genes_mitochondrial.tsv'
    log:
        "logs/genes/genes_mitochondrial.log"
    script:
        "../../scripts/mitochondrial_genes.py"

rule curated_markers_counts_spatial:
    input:
        tsv='data/curated_markers.tsv'
    output:
        dir=directory('figures/spatial/markers_counts/')
    params:
        samples=config["samples"],
    log:
        "logs/curated_markers_counts_spatial.log"
    conda:
        "../../envs/scanpy-env.yaml"
    script:
        "../../scripts/curated_markers_counts_spatial.py"
