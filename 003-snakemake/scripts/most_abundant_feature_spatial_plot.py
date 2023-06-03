# Based on <https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html?highlight=visium#1.-Loading-Visium-data>

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
for i in sample_names:
    slides.append(read_and_qc(i))

# Combine anndata objects together
adata = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_names,
    index_unique=None
)
#######################

##
# Fetch list of genes to visualise
##

input_genes_table = pd.read_table(snakemake.input["tsv"], sep="\t")
# remove ribosomal and mitochondrial to get to the interesting stuff
filtered_genes_table = input_genes_table[[not x.startswith(("RPS", "RPL", "MT-")) for x in input_genes_table['SYMBOL']]]
# peek at first 8 genes (2 rows of 4 plots)
gene_symbols = filtered_genes_table['SYMBOL'].tolist()[:8]

for sample_name in sample_names:
    slide = select_slide(adata, sample_name)
    with mpl.rc_context({'figure.figsize': [6,7],
                        'axes.facecolor': 'black'}):
        sc.pl.spatial(slide,
                    color=gene_symbols, img_key=None, size=1,
                    vmin=0, cmap='magma', vmax='p90.0',
                    gene_symbols='SYMBOL', save=f"-{sample_name}.png"
                    )
    os.rename(f"show-{sample_name}.png", f"figures/spatial/most_detected_most_abundant_features/{sample_name}.png")
