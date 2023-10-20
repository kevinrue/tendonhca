rule seurat_transfer:
    input:
        "results/spaceranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        predictions_png="figures/seurat_transfer/predictions/{sample}.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_transfer/{sample}.log",
    threads: 1
    resources:
        mem_mb=20 * 1024,
        runtime="20m",
    script:
        "../../scripts/seurat_transfer_data.R"