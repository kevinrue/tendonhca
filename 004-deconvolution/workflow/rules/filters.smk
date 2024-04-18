rule cellranger_filter_barcodes:
    input:
        h5="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        barcodes="results/cellranger_filter_barcodes/{sample}.tsv.gz",
    params:
        umi_cutoff_min=lambda wildcards, output: get_umi_min_cutoff(wildcards),
    conda:
        "../../envs/r-env.yaml"
    log:
        "results/cellranger_filter_barcodes/{sample}.log",
    script:
        "../../scripts/cellranger_filter_barcodes.R"
