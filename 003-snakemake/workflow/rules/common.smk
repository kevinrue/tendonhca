import yaml
import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_final_output():
    final_output = ["results/spaceranger_stats/runtime.tsv",
                    "results/spaceranger_stats/total_counts.tsv"]
    final_output.append(expand(
        "results/jupyter_notebooks/spaceranger_qc/{sample}.html",
        sample=samples.index.tolist(),
    ))
    return final_output

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}

def get_image(wildcards):
    u = samples.loc[wildcards.sample]
    return {"image": u["image"]}

def get_slide(wildcards):
    u = samples.loc[wildcards.sample]
    return u["slide"]

def get_area(wildcards):
    u = samples.loc[wildcards.sample]
    return u["area"]