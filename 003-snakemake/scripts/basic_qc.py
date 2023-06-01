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
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]
    # TODO: identify ribosomal genes
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    # from scipy.sparse import csr_matrix
    # adata.X = adata.X.toarray()
    # adata.X = csr_matrix(adata.X)
    # adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    
    return adata

def select_slide(adata, s, s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """
    
    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    
    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}
    
    return slide

#######################

# Read the data into anndata objects
slides = []
for i in [sample_name]:
    slides.append(read_and_qc(i))

# Combine anndata objects together
adata = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=[sample_name],
    index_unique=None
)
#######################

# mitochondria-encoded (MT) genes should be removed for spatial mapping
# we do not do this for this basic quality control, to include everything
# adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
# adata = adata[:, ~adata.var['mt'].values]

##
# plot quality control metrics for each sample
##

fig, axs = plt.subplots(2, 4, figsize=(15, 8))
for i, s in enumerate(adata.obs['sample'].unique()):
    #fig.suptitle('Covariates for filtering')
    slide = select_slide(adata, s)
    
    sns.distplot(slide.obs['total_counts'],
                 kde=False, ax = axs[0, 0])
    axs[0, 0].set_xlim(0, adata.obs['total_counts'].max())
    axs[0, 0].set_xlabel(f'total_counts | {s}')
    
    x_max = np.quantile(slide.obs['total_counts'], .9)
    sns.distplot(slide.obs['total_counts']\
                 [slide.obs['total_counts']<x_max],
                 kde=False, bins=40, ax = axs[0, 1])
    axs[0, 1].set_xlim(0, x_max)
    axs[0, 1].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts'],
                 kde=False, bins=60, ax = axs[0, 2])
    axs[0, 2].set_xlim(0, adata.obs['n_genes_by_counts'].max())
    axs[0, 2].set_xlabel(f'n_genes_by_counts | {s}')
    
    x_max = np.quantile(slide.obs['n_genes_by_counts'], .9)
    sns.distplot(slide.obs['n_genes_by_counts']\
                 [slide.obs['n_genes_by_counts']<x_max],
                 kde=False, bins=60, ax = axs[0, 3])
    axs[0, 3].set_xlim(0, x_max)
    axs[0, 3].set_xlabel(f'n_genes_by_counts | {s}')
    
    sns.distplot(slide.obs['pct_counts_mt'],
                 kde=False, ax = axs[1, 0])
    axs[1, 0].set_xlim(0, adata.obs['pct_counts_mt'].max())
    axs[1, 0].set_xlabel(f'pct_counts_mt | {s}')
    
    x_max = np.quantile(slide.obs['pct_counts_mt'], .9)
    sns.distplot(slide.obs['pct_counts_mt']\
                 [slide.obs['pct_counts_mt']<x_max],
                 kde=False, bins=40, ax = axs[1, 1])
    axs[1, 1].set_xlim(0, adata.obs['pct_counts_mt'].quantile(0.9))
    axs[1, 1].set_xlabel(f'pct_counts_mt | {s}')

plt.savefig(snakemake.output['histogram'])

##
# plot quality control metrics in spatial coordinates
##

slide = select_slide(adata, sample_name)

with mpl.rc_context({'figure.figsize': [6,7],
                     'axes.facecolor': 'white'}):
    print(sc.settings.figdir)
    fig = sc.pl.spatial(slide, img_key = "hires", cmap='magma', ncols=2,
                  library_id=list(slide.uns['spatial'].keys())[0],
                  color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], size=1,
                  vmin=0, vmax='p90.0',
                  gene_symbols='SYMBOL', show=False, return_fig=True,
                  save=f"-sc_pl_spatial-{sample_name}.png")
    os.rename(f"show-sc_pl_spatial-{sample_name}.png", snakemake.output["spatial"])

##
# identify most abundant genes
##

adata_feature_mean = np.array(adata.X.mean(axis=0)).flatten()
d = {'SYMBOL': adata.var['SYMBOL'], 'mean': adata_feature_mean}
adata_feature_mean = pd.DataFrame(d, index=adata.var_names).sort_values(by=['mean'], ascending=False)
adata_feature_mean.head(100).to_csv(snakemake.output["features_mean_top_100"], sep="\t")
