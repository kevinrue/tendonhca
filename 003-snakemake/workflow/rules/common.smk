import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_final_output():
    final_output = expand(
        "results/spaceranger_count/{sample}/outs/filtered_feature_bc_matrix.h5",
        sample=samples["sample_name"],
    )
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