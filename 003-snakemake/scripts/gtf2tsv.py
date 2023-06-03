##
# imports
##

import re
import pandas as pd

##
# main
##

genes_dict = {}
with open(snakemake.input['gtf'], 'rt') as stream_in:
    for line in stream_in:
        if line.startswith('#'):
            continue
        else:
            line_data = line.strip().split('\t')
            if line_data[2] == 'gene':
                gene_name = re.search(r'gene_name \"(.+?)\";', line_data[8]).group(1)
                gene_id = re.search(r'gene_id \"(.+?)\";', line_data[8]).group(1)
                genes_dict[gene_id] = {'gene_name': gene_name, 'seqname': line_data[0]}

##
# outputs
##

genes_df = pd.DataFrame.from_dict(genes_dict, orient='index')
genes_df.to_csv(snakemake.output['tsv'], sep='\t', index_label = 'gene_id')
