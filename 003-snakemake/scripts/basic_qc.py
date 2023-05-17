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
print(sc.settings.figdir)

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
    
    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    
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
adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
adata = adata[:, ~adata.var['mt'].values]

# PLOT QC FOR EACH SAMPLE
fig, axs = plt.subplots(len(slides), 4, figsize=(15, 4))
for i, s in enumerate(adata.obs['sample'].unique()):
    #fig.suptitle('Covariates for filtering')
    slide = select_slide(adata, s)
    sns.distplot(slide.obs['total_counts'],
                 kde=False, ax = axs[0])
    axs[0].set_xlim(0, adata.obs['total_counts'].max())
    axs[0].set_xlabel(f'total_counts | {s}')
    x_max = np.quantile(slide.obs['total_counts'], .9)
    sns.distplot(slide.obs['total_counts']\
                 [slide.obs['total_counts']<x_max],
                 kde=False, bins=40, ax = axs[1])
    axs[1].set_xlim(0, x_max)
    axs[1].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts'],
                 kde=False, bins=60, ax = axs[2])
    axs[2].set_xlim(0, adata.obs['n_genes_by_counts'].max())
    axs[2].set_xlabel(f'n_genes_by_counts | {s}')
    x_max = np.quantile(slide.obs['n_genes_by_counts'], .9)
    sns.distplot(slide.obs['n_genes_by_counts']\
                 [slide.obs['n_genes_by_counts']<x_max],
                 kde=False, bins=60, ax = axs[3])
    axs[3].set_xlim(0, x_max)
    axs[3].set_xlabel(f'n_genes_by_counts | {s}')

plt.savefig(snakemake.output['total_counts_n_genes_by_counts'])

## plot quality control metrics in spatial coordinates

slide = select_slide(adata, sample_name)

with mpl.rc_context({'figure.figsize': [6,7],
                     'axes.facecolor': 'white'}):
    print(sc.settings.figdir)
    fig = sc.pl.spatial(slide, img_key = "hires", cmap='magma',
                  library_id=list(slide.uns['spatial'].keys())[0],
                  color=['total_counts', 'n_genes_by_counts'], size=1,
                  vmin=0, vmax='p90.0',
                  gene_symbols='SYMBOL', show=False, return_fig=True,
                  save=f"-sc_pl_spatial-{sample_name}.png")
    os.rename(f"show-sc_pl_spatial-{sample_name}.png", snakemake.output["total_counts_n_genes_by_counts_spatial"])
