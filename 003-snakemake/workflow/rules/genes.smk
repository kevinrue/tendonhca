rule gene_table:
    input:
        gtf=config['transcriptome'] + 'genes/genes.gtf'
    output:
        tsv='annotations/genes.tsv'
    script:
        "../../scripts/gtf2tsv.py"
