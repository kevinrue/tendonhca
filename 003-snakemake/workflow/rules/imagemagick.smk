rule spatial_basic_qc:
    input:
        unpack(get_image),
        qc="results/qc/total_counts_n_genes_by_counts_spatial/{sample}.png",
    output:
        png="results/qc/total_counts_n_genes_by_counts_spatial/{sample}_slide.png",
    shell:
        "montage "
        "{input.image} "
        "{input.qc} "
        "-geometry 1000x600\> "
        "{output.png}"
