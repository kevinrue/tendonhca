import scanpy as sc
import pandas as pd

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

path = "results/spaceranger_count/"

samples = (
    pd.read_csv(snakemake.config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

sample_names = samples.index.tolist()

with open(snakemake.output[0], "wt") as filestream_out:
    filestream_out.write("\t".join(["sample_name", "total_counts"]) + "\n")
    for sample_name in sample_names:
        adata = sc.read_visium(path + str(sample_name) + '/outs',
                               count_file='filtered_feature_bc_matrix.h5', load_images=True)
        # rename observations and variables
        adata.obs['sample'] = sample_name
        adata.var['SYMBOL'] = adata.var_names
        adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
        adata.var_names = adata.var['ENSEMBL']
        adata.var.drop(columns='ENSEMBL', inplace=True)
        # compute qc metrics
        adata.X = adata.X.toarray()
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        total_counts = int(sum(adata.obs['total_counts']))
        filestream_out.write(f"{sample_name}\t{total_counts}\n")
