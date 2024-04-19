rule vireo_qc:
    input:
        cellranger_filtered="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        vireo_dir="results/vireo/{sample}"
    output:
        donor_umi="results/vireo_qc/{sample}/donor_umi.pdf",
        donor_n_vars="results/vireo_qc/{sample}/donor_n_vars.pdf",
        prob_donor_doublet="results/vireo_qc/{sample}/prob_donor_doublet.pdf",
    params:
        sample="{sample}"
    conda:
        "../../envs/r-env.yaml"
    log:
        "results/vireo_qc/{sample}.log",
    script:
        "../../scripts/vireo_qc.R"

rule vireo_filtered_barcodes_qc:
    input:
        cellranger_filtered="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        vireo_dir="results/vireo_filtered_barcodes/{sample}"
    output:
        donor_umi="results/vireo_filtered_barcodes_qc/{sample}/donor_umi.pdf",
        donor_n_vars="results/vireo_filtered_barcodes_qc/{sample}/donor_n_vars.pdf",
        prob_donor_doublet="results/vireo_filtered_barcodes_qc/{sample}/prob_donor_doublet.pdf",
    params:
        sample="{sample}"
    conda:
        "../../envs/r-env.yaml"
    log:
        "results/vireo_filtered_barcodes_qc/{sample}.log",
    script:
        "../../scripts/vireo_qc.R"

rule vireo_qc_sample_cowplot:
    input:
        donor_umi="results/vireo_qc/{sample}/donor_umi.pdf",
        donor_n_vars="results/vireo_qc/{sample}/donor_n_vars.pdf",
        prob_donor_doublet="results/vireo_qc/{sample}/prob_donor_doublet.pdf"
    output:
        montage="results/vireo_qc_cowplot/{sample}.pdf"
    conda:
        "../../envs/r-env.yaml"
    script:
        "../../scripts/vireo_qc_sample_cowplot.R"

rule vireo_filtered_barcodes_qc_sample_cowplot:
    input:
        donor_umi="results/vireo_filtered_barcodes_qc/{sample}/donor_umi.pdf",
        donor_n_vars="results/vireo_filtered_barcodes_qc/{sample}/donor_n_vars.pdf",
        prob_donor_doublet="results/vireo_filtered_barcodes_qc/{sample}/prob_donor_doublet.pdf"
    output:
        montage="results/vireo_filtered_barcodes_qc_cowplot/{sample}.pdf"
    conda:
        "../../envs/r-env.yaml"
    script:
        "../../scripts/vireo_qc_sample_cowplot.R"
