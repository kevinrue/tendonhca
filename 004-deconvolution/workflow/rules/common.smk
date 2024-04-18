import yaml
import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_final_output():
    final_output = []
    final_output.append(expand(
        "results/vireo_qc_cowplot/{sample}.pdf",
        sample=samples.index.tolist(),
    ))
    final_output.append(expand(
        "results/cellranger_filter_barcodes/{sample}.txt",
        sample=samples.index.tolist(),
    ))
    return final_output

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}

def get_n_donors(wildcards):
    u = samples.loc[wildcards.sample]
    return {"n_donors": u["n_donors"]}

def get_umi_min_cutoff(wildcards):
    u = samples.loc[wildcards.sample]
    return {"umi_cutoff_min": u["umi_cutoff_min"]}
