rule spatial_basic_qc:
    input:
        unpack(get_image),
        qc="figures/initial_qc/spatial/{sample}/metrics.png",
    output:
        png="figures/initial_qc/spatial/{sample}/metrics_slide.png",
    shell:
        "montage "
        "{input.image} "
        "{input.qc} "
        "-geometry 600x600\> "
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
        "-geometry 600x600\> "
        "{output.png}"

rule spatial_clusters:
    input:
        unpack(get_image),
        feature="figures/filtered_genes/{sample}/dimred/spatial_clusters.png",
    output:
        png="figures/filtered_genes/{sample}/dimred/spatial_clusters_image.png",
    shell:
        "montage "
        "{input.image} "
        "{input.feature} "
        "-geometry 600x600\> "
        "{output.png}"
