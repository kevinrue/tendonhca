import pandas as pd

genes_df = pd.read_csv(snakemake.input['tsv'], sep='\t')

ribosomal_df = genes_df[[gene_name.startswith(('RPS', 'RPL')) for gene_name in genes_df["gene_name"]]]\
    [['gene_id', 'gene_name']]

ribosomal_df.to_csv(snakemake.output['tsv'], sep='\t', index=False)
