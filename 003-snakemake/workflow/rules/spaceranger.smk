# based on https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/e103c1cc78feba97cc3cebe8d7f2a51c8958ab96/workflow/rules/align.smk
rule spaceranger:
    input:
        directory(unpack(get_fastqs)),
        unpack(get_image),
    output:
        filtered_feature_bc_matrix="results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    log:
        "logs/spaceranger_count/{sample}.log",
    params:
        slide=lambda wildcards, output: get_slide(wildcards),
        area=lambda wildcards, output: get_area(wildcards),
        transcriptome=config["transcriptome"],
        localcores=config["localcores"],
        localmem=config["localmem"],
    threads: 6
    resources:
        mem_mb=config["localmem"] * 1000,
        runtime=config["runtime"],
    envmodules:
        "spaceranger/1.3.1",
    shell:
        "spaceranger count --id={wildcards.sample} "
        "--transcriptome={params.transcriptome} "
        "--fastqs={input.fastqs} "
        "--sample={wildcards.sample} "
        "--image={input.image} "
        "--slide={params.slide} "
        "--area={params.area} "
        "--localcores={params.localcores} "
        "--localmem={params.localmem} && "
        "rm -rf results/spaceranger_count/{wildcards.sample} && "
        "mv {wildcards.sample} results/spaceranger_count/"
