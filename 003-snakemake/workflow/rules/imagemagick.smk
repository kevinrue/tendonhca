rule spatial_basic_qc:
    input:
        unpack(get_image),
        qc="figures/basic_qc/spatial/metrics/{sample}.png",
    output:
        png="figures/basic_qc/spatial/slide/{sample}.png",
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


