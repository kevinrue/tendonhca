rule seurat_transfer:
    input:
        "results/spaceranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        predictions_png="figures/seurat/transfer/predictions/{sample}.png",
        top_prediction_png="figures/seurat/transfer/top_prediction/{sample}.png",
        probabilities_diff_png="figures/seurat/transfer/predictions/histogram_probability_diff_{sample}.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_transfer/{sample}.log",
    threads: 2
    resources:
        mem_mb=20 * 1024,
        runtime="20m",
    script:
        "../../scripts/seurat_transfer_data.R"

rule seurat_qc_sample:
    input:
        "results/spaceranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        histograms_png="figures/seurat_qc/{sample}_histogram.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_qc_sample/{sample}_histogram.log",
    threads: 2
    resources:
        mem_mb=4 * 1024,
        runtime="5m",
    script:
        "../../scripts/seurat_qc_sample.R"

rule seurat_transfer_montage:
    input:
        lowres_png="results/spaceranger_count/{sample}/outs/spatial/tissue_lowres_image.png",
        predictions_png="figures/seurat/transfer/predictions/{sample}.png",
        probabilities_diff_png="figures/seurat/transfer/predictions/histogram_probability_diff_{sample}.png",
        qcmetrics_png="figures/seurat_qc/{sample}_histogram.png",
        slide_top_prediction_png="figures/seurat/transfer/top_prediction/{sample}.png",
    output:
        montage_png="figures/seurat/predictions_montage/{sample}.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_transfer_montage/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1024,
        runtime="5m",
    script:
        "../../scripts/cowplot_seurat_predictions.R"

# different reference for each sample
rule seurat_transfer_reference:
    input:
        "results/spaceranger_count/{sample}/outs/raw_feature_bc_matrix.h5",
    output:
        predictions_png="figures/seurat/transfer_reference/predictions/{sample}.png",
        top_prediction_png="figures/seurat/transfer_reference/top_prediction/{sample}.png",
        probabilities_diff_png="figures/seurat/transfer_reference/predictions/histogram_probability_diff_{sample}.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_transfer_reference/{sample}.log",
    params:
        samples=config["samples"],
    threads: 2
    resources:
        mem_mb=80 * 1024,
        runtime="20m",
    script:
        "../../scripts/seurat_transfer_data_reference.R"

rule seurat_transfer_montage_reference:
    input:
        lowres_png="results/spaceranger_count/{sample}/outs/spatial/tissue_lowres_image.png",
        predictions_png="figures/seurat/transfer_reference/predictions/{sample}.png",
        probabilities_diff_png="figures/seurat/transfer_reference/predictions/histogram_probability_diff_{sample}.png",
        qcmetrics_png="figures/seurat_qc/{sample}_histogram.png",
        slide_top_prediction_png="figures/seurat/transfer_reference/top_prediction/{sample}.png",
    output:
        montage_png="figures/seurat/transfer_reference/montage/{sample}.png",
    conda:
        "../../envs/r-env.yaml"
    log:
        "logs/seurat_transfer_montage_reference/{sample}.log",
    threads: 1
    resources:
        mem_mb=2 * 1024,
        runtime="5m",
    script:
        "../../scripts/cowplot_seurat_predictions_reference.R"
