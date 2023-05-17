rule basic_qc:
    input:
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
    output:
        "results/jupyter_notebooks/spaceranger_qc/{sample}.html",
    conda:
        "../../envs/scanpy-env.yaml"
    log:
        # optional path to the processed notebook
        notebook="logs/notebooks/basic_qc/{sample}.ipynb"
    threads: 1
    resources:
        mem_mb=2 * 1000,
        runtime="5m",
    notebook:
        "../../notebooks/basic_qc.py.ipynb"

    