rule spatial_basic_qc:
    input:
        unpack(get_image),
        qc="figures/qc/total_counts_n_genes_by_counts_spatial/{sample}.png",
    output:
        png="figures/qc/total_counts_n_genes_by_counts_spatial/{sample}_slide.png",
    shell:
        "montage "
        "{input.image} "
        "{input.qc} "
        "-geometry 1000x600\> "
        "{output.png}"

rule spatial_most_detected:
    input:
        unpack(get_image),
        feature="figures/spatial/most_detected_most_abundant_features/{sample}.png",
    output:
        png="figures/spatial/most_detected_most_abundant_features/{sample}_image.png",
    shell:
        "montage "
        "{input.image} "
        "{input.feature} "
        "-geometry 1000x600\> "
        "{output.png}"


