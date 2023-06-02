import pandas as pd

genes_df = pd.read_csv(snakemake.input['tsv'], sep='\t')

mt_df = genes_df[genes_df["seqname"] == "chrM"][['gene_id', 'gene_name']]

mt_df.to_csv(snakemake.output['tsv'], sep='\t', index=False)
