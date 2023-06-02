import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

# change default output directory for figures
sc.settings.figdir = './'

sample_name = snakemake.wildcards['sample']

mitochondrial_genes = pd.read_csv(snakemake.input['mitochondrial'], sep="\t")
ribosomal_genes = pd.read_csv(snakemake.input['ribosomal'], sep="\t")

def read_and_qc(sample_name):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.
    
    :param sample_name: Name of the sample
    """
    
    adata = sc.read_visium("results/spaceranger_count/" + str(sample_name) + '/outs',
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)
    
    # some filtered matrices have spots with no counts (!), so we need to remove them
    sc.pp.filter_cells(adata, min_counts=1)
    # identify mitochondria-encoded genes
    adata.var['mt'] = [symbol in set(mitochondrial_genes['gene_name']) for symbol in adata.var['SYMBOL']]
    adata = adata[:, ~adata.var['mt'].values]
    # identify ribosomal genes
    adata.var['ribosomal'] = [symbol in set(ribosomal_genes['gene_name']) for symbol in adata.var['SYMBOL']]
    adata = adata[:, ~adata.var['ribosomal'].values]
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    
    s_keys = list(adata.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[sample_name in k for k in s_keys]][0]
    
    adata.uns['spatial'] = {s_spatial: adata.uns['spatial'][s_spatial]}
    
    return adata

slide = read_and_qc(sample_name)

##
# plot quality control metrics for each sample
##

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
#fig.suptitle('Covariates for filtering')
    
sns.distplot(slide.obs['total_counts'],
                kde=False, ax = axs[0])
axs[0].set_xlim(0, slide.obs['total_counts'].max())
axs[0].set_xlabel(f'total_counts | {sample_name}')

x_max = slide.obs['total_counts'].quantile(0.9)
sns.distplot(slide.obs['total_counts']\
                [slide.obs['total_counts']<x_max],
                kde=False, bins=40, ax = axs[1])
axs[1].set_xlim(0, x_max)
axs[1].set_xlabel(f'total_counts | {sample_name}')

sns.distplot(slide.obs['n_genes_by_counts'],
                kde=False, bins=60, ax = axs[2])
axs[2].set_xlim(0, slide.obs['n_genes_by_counts'].max())
axs[2].set_xlabel(f'n_genes_by_counts | {sample_name}')

x_max = slide.obs['n_genes_by_counts'].quantile(0.9)
sns.distplot(slide.obs['n_genes_by_counts']\
                [slide.obs['n_genes_by_counts']<x_max],
                kde=False, bins=60, ax = axs[3])
axs[3].set_xlim(0, x_max)
axs[3].set_xlabel(f'n_genes_by_counts | {sample_name}')

plt.savefig(snakemake.output['histogram'])

##
# plot quality control metrics in spatial coordinates
##

with mpl.rc_context({'figure.figsize': [6,3.5],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.spatial(slide, img_key = "hires", cmap='magma', ncols=2,
                  library_id=list(slide.uns['spatial'].keys())[0],
                  color=['total_counts', 'n_genes_by_counts'], size=1,
                  vmin=0, vmax='p90.0',
                  gene_symbols='SYMBOL', show=False, return_fig=True,
                  save=f"-sc_pl_spatial-{sample_name}.png")
    os.rename(f"show-sc_pl_spatial-{sample_name}.png", snakemake.output["spatial"])

##
# identify most abundant genes
##

adata_feature_mean = np.array(slide.X.mean(axis=0)).flatten()
d = {'SYMBOL': slide.var['SYMBOL'], 'mean': adata_feature_mean}
adata_feature_mean = pd.DataFrame(d, index=slide.var_names).sort_values(by=['mean'], ascending=False)
adata_feature_mean.head(100).to_csv(snakemake.output["features_mean_top_100"], sep="\t")
