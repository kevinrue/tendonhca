rule cellranger_filter_barcodes:
    input:
        h5="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        barcodes="results/cellranger_filter_barcodes/{sample}.txt",
    params:
        slide=lambda wildcards, output: get_umi_min_cutoff(wildcards),
    log:
        "results/cellranger_filter_barcodes/{sample}.log",
    script:
        "../../scripts/cellranger_filter_barcodes.R"
