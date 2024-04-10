import yaml
import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_final_output():
    final_output = []
    # initial quality control
    final_output.append(expand(
        "results/cellranger_count/{sample}/outs/possorted_genome_bam.bam",
        sample=samples.index.tolist(),
    ))
    return final_output

def get_fastqs(wildcards):
    u = samples.loc[wildcards.sample]
    return {"fastqs": u["fastqs"]}
