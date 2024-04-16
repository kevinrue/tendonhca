rule cellsnp_lite:
    input:
        dir="results/cellsnp-lite/{sample}",
    output:
        summary="results/vireo/{sample}/summary.tsv",
    log:
        "results/cellranger_count/{sample}.log",
    params:
        slide=lambda wildcards, output: get_n_donors(wildcards),
    shell:
        "vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR"
        "-c {input.dir}"
        "-N {params.n_donor}"
        "-o {output.summary}"

