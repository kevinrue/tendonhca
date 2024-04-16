rule vireo_qc:
    input:
        cellranger_filtered="results/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        vireo_donor_ids="results/vireo/{sample}/donor_ids.tsv"
    output:
        pdf="results/vireo_qc/{sample}/donor_umi.pdf",
    conda:
        "../../envs/r-env.yaml"
    log:
        "results/vireo_qc/{sample}.log",
    script:
        "../../scripts/vireo_qc.R"

