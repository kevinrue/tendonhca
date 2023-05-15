# based on https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/e103c1cc78feba97cc3cebe8d7f2a51c8958ab96/workflow/rules/align.smk
rule spaceranger:
    input:
        unpack(get_fastqs),
    output:
        filtered_feature_bc_matrix="results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    log:
        "logs/spaceranger_count/{sample}.log",
    threads: 6
    resources:
        mem_mb: 25000
    shell:
        "module load spaceranger/1.3.1 &&"
        "spaceranger count --id={wildcard.sample}"
        "--transcriptome=/project/tendonhca/shared/spatial/analysis/refdata-gex-GRCh38-2020-A/"
        "--fastqs={fastqs}"
        "--sample={wildcard.sample}"
        "--image=/project/tendonhca/shared/spatial/analysis/achilles/images/OMB1556_Ach_Enth.jpeg"
        "--slide=V12J03-133"
        "--area=C1"
        "--localcores=6"
        "--localmem=25"
