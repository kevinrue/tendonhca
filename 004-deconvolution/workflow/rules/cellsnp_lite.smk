rule cellsnp_lite:
    input:
        bam="results/cellsnp_lite/{sample}/outs/possorted_genome_bam.bam",
        barcode="results/cellsnp_lite/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        vcf="results/cellsnp-lite/{sample}/cellSNP.base.vcf.gz",
    log:
        "results/cellranger_count/{sample}.log",
    shell:
        "cellsnp-lite "
        "-s {input.bam}"
        "-b {input.barcode}"
        "-O results/cellsnp-lite/{wildcards.sample} "
        "-p 10"
        "--minMAF 0.1"
        "--minCOUNT 100"
        "--gzip"
