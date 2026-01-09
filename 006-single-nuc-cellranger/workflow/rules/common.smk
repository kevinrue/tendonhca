import yaml
import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}
