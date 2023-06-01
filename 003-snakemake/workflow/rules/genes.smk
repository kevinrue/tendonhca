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

rule curated_markers_counts_spatial_full:
    input:
        tsv='data/curated_markers.tsv'
    output:
        dir=directory('figures/spatial/curated_celltype_markers_counts_full'),
    params:
        samples=config["samples"],
    log:
        "logs/curated_markers_counts_spatial_full.log"
    conda:
        "../../envs/scanpy-env.yaml"
    script:
        "../../scripts/curated_markers_counts_spatial_full.py"

rule curated_markers_counts_spatial_quantile:
    input:
        tsv='data/curated_markers.tsv'
    output:
        dir=directory('figures/spatial/curated_celltype_markers_counts_quantile'),
    params:
        samples=config["samples"],
    log:
        "logs/curated_markers_counts_spatial_quantile.log"
    conda:
        "../../envs/scanpy-env.yaml"
    script:
        "../../scripts/curated_markers_counts_spatial_quantile.py"

rule curated_markers_counts_spatial_log1p:
    input:
        tsv='data/curated_markers.tsv'
    output:
        dir=directory('figures/spatial/curated_celltype_markers_counts_log1p'),
    params:
        samples=config["samples"],
    log:
        "logs/curated_markers_counts_spatial_log1p.log"
    conda:
        "../../envs/scanpy-env.yaml"
    script:
        "../../scripts/curated_markers_counts_spatial_log1p.py"

