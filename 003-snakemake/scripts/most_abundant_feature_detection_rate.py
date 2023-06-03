import glob
import os
import pandas as pd

input_files = glob.glob("results/qc/features_mean_top_100/*.tsv")

input_tables = []
for input_file in input_files:
    input_table = pd.read_table(input_file, sep="\t")
    sample_name = os.path.splitext(os.path.basename(input_file))[0]
    input_table['sample_name'] = sample_name
    input_tables.append(input_table)

concat_input_tables = pd.concat(input_tables)
tally = concat_input_tables.groupby(['ENSEMBL', 'SYMBOL']).size().sort_values(ascending=False) / len(input_files)

df_tally = pd.DataFrame({'tally': tally})
df_tally.to_csv(snakemake.output["tsv"], sep="\t")

##
# extras
##

# identify top genes that are not ribosomal
# df_tally[[not x.startswith(("RPS", "RPL")) for x in df_tally['SYMBOL']]]
