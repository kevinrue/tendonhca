rule cellsnp_lite:
    input:
        bam="results/cellsnp_lite/{sample}/outs/possorted_genome_bam.bam",
        barcode="results/cellsnp_lite/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    output:
        dir=directory("results/cellsnp-lite/{sample}"),
    log:
        "results/cellsnp-lite/{sample}.log",
    conda:
        "../../envs/cellsnp_lite.yaml",
    threads: 10
    resources:
        mem_mb=16 * 1024,
        runtime="23h",
    shell:
        "cellsnp-lite "
        "-s {input.bam} "
        "-b {input.barcode} "
        "-O {output.dir} "
        "-p 10 "
        "--minMAF 0.1 "
        "--minCOUNT 100 "
        "--gzip "

rule cellsnp_lite_filtered_barcodes:
    input:
        bam="results/cellranger_count/{sample}/outs/possorted_genome_bam.bam",
        barcode="results/cellranger_filter_barcodes/{sample}.tsv.gz",
    output:
        dir=directory("results/cellsnp-lite_filtered_barcodes/{sample}"),
    log:
        "results/cellsnp-lite_filtered_barcodes/{sample}.log",
    conda:
        "../../envs/cellsnp_lite.yaml",
    threads: 10
    resources:
        mem_mb=16 * 1024,
        runtime="23h",
    shell:
        "cellsnp-lite "
        "-s {input.bam} "
        "-b {input.barcode} "
        "-O {output.dir} "
        "-p 10 "
        "--minMAF 0.1 "
        "--minCOUNT 100 "
        "--gzip "
