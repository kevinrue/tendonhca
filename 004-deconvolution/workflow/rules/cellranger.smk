rule cellranger_count:
    input:
        unpack(get_fastqs),
    output:
        bam="results/cellranger_count/{sample}/outs/possorted_genome_bam.bam",
    log:
        "results/cellranger_count/{sample}.log",
    params:
        transcriptome=config["transcriptome"],
        localcores=config["localcores"],
        localmem=config["localmem"],
    threads: config["localcores"]
    resources:
        mem_mb=config["localmem"] * 1024,
        runtime=config["runtime"],
    envmodules:
        "cellranger/7.2.0",
    shell:
        "cellranger count "
        "--id={wildcards.sample} "
        "--transcriptome={params.transcriptome} "
        "--fastqs={input.fastqs} "
        "--sample={wildcards.sample} "
        "--create-bam=true "
        "--localcores={params.localcores} "
        "--localmem={params.localmem} && "
        "rm -rf results/cellranger_count/{wildcards.sample} && "
        "mv {wildcards.sample} results/cellranger_count/"
