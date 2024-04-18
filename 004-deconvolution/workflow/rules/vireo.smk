rule vireo:
    input:
        dir="results/cellsnp-lite/{sample}",
    output:
        summary="results/vireo/{sample}/summary.tsv",
    log:
        "results/vireo/{sample}.log",
    params:
        n_donors=lambda wildcards, output: get_n_donors(wildcards),
    conda:
        "../../envs/vireo.yaml"
    shell:
        "vireo "
        "-c {input.dir}"
        "-N {params.n_donors}"
        "-o {output.summary}"

rule vireo_filtered_barcodes:
    input:
        dir="results/cellsnp-lite_filtered_barcodes/{sample}",
    output:
        summary="results/vireo_filtered_barcodes/{sample}/summary.tsv",
    log:
        "results/vireo_filtered_barcodes/{sample}.log",
    params:
        n_donors=lambda wildcards, output: get_n_donors(wildcards),
    conda:
        "../../envs/vireo.yaml"
    shell:
        "vireo "
        "-c {input.dir} "
        "-N {params.n_donors} "
        "-o {output.summary} "
