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

##
# import information about samples
##

samples = (
    pd.read_csv(snakemake.config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

sample_names = samples.index.tolist()

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
    
    # some filtered matrices have spots with no counts (!), so we need to remove them before computing QC metrics
    sc.pp.filter_cells(adata, min_counts=samples['min_counts'][sample_name])
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
# Process data
##

sc.pp.normalize_total(slide) # target_sum = median of counts per cell
sc.pp.log1p(slide)
sc.pp.highly_variable_genes(slide, flavor="seurat", n_top_genes=2000) # TODO: set n_top_genes per sample
sc.pp.pca(slide, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(slide, n_neighbors=15, n_pcs=40)
sc.tl.umap(slide)
sc.tl.leiden(slide, key_added="clusters")

##
# Produce plots
##

with mpl.rc_context({'figure.figsize': [6,6],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.highly_variable_genes(slide)
    plt.savefig(snakemake.output['highly_variable_genes'])

with mpl.rc_context({'figure.figsize': [6,6],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.pca_variance_ratio(slide, log=True)
    plt.savefig(snakemake.output['pca_variance_ratio'])

with mpl.rc_context({'figure.figsize': [6,6],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.pca(slide, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
    plt.savefig(snakemake.output['pca'])

with mpl.rc_context({'figure.figsize': [6,6],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.umap(slide, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
    plt.savefig(snakemake.output['umap'])

with mpl.rc_context({'figure.figsize': [6,6],
                     'axes.facecolor': 'white'}):
    fig = sc.pl.spatial(slide, img_key="hires", color="clusters", size=1.5)
    plt.savefig(snakemake.output['spatial_clusters'])
