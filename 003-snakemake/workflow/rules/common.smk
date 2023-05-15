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

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}
