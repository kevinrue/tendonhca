rule cellsnp_lite:
    input:
        bam="results/cellsnp_lite/{sample}/outs/possorted_genome_bam.bam",
        barcode="results/cellsnp_lite/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        dir=directory("results/cellsnp-lite/{sample}"),
    log:
        "results/cellsnp_lite/{sample}.log",
    shell:
        "cellsnp-lite "
        "-s {input.bam}"
        "-b {input.barcode}"
        "-O results/cellsnp-lite/{output.dir} "
        "-p 10"
        "--minMAF 0.1"
        "--minCOUNT 100"
        "--gzip"

rule cellsnp_lite_filtered_barcodes:
    input:
        bam="results/cellsnp_lite_filtered_barcodes/{sample}/outs/possorted_genome_bam.bam",
        barcode="results/cellsnp_lite_filtered_barcodes/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        dir=directory("results/cellsnp-lite/{sample}"),
    log:
        "results/cellsnp_lite_filtered_barcodes/{sample}.log",
    shell:
        "cellsnp-lite "
        "-s {input.bam}"
        "-b {input.barcode}"
        "-O results/cellsnp-lite/{output.dir} "
        "-p 10"
        "--minMAF 0.1"
        "--minCOUNT 100"
        "--gzip"
