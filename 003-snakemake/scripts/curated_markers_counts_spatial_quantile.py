##
# imports
##

import os
import re

import IPython

import pandas as pd

import scanpy as sc

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

# change default output directory for figures
sc.settings.figdir = './'

if os.path.exists(snakemake.output["dir"]):
    os.rmdir(snakemake.output["dir"])

##
# import information about samples
##

samples = (
    pd.read_csv(snakemake.config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

sample_names = samples.index.tolist()

##
# define functions
##

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

# TODO: process each sample separately to minimise memory footprint

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
# Load curated markers
##

curated_markers = pd.read_csv(snakemake.input["tsv"], sep="\t")
curated_markers.head()

##
# Plot markers next to total counts
##

# NOTE: 95% quantile used to cap color scale (i.e., trim outlier spots with high expression)

for sample_name in sample_names:
    # fetch data for the sample
    slide = select_slide(adata, sample_name)
    sc.pp.filter_cells(slide, min_counts=100)
    for cell_type in curated_markers['cell_type'].unique():
        # fetch markers for the cell type
        markers_symbols = curated_markers['gene_symbol'][curated_markers['cell_type'] == cell_type].tolist()
        # TODO: some manually curated markers are not present in the data
        markers_symbols_available = list(set(markers_symbols).intersection(slide.var['SYMBOL']))
        # remove spaces prior to generating file names
        cell_type = cell_type.replace(" ", "_")
        # save png file to current directory (temporarily)
        fig = sc.pl.spatial(
            slide, img_key = "hires", cmap='magma',
            library_id=list(slide.uns['spatial'].keys())[0],
            color=['total_counts'] + markers_symbols_available, size=1,
            vmin=0, vmax='p95.0',
            gene_symbols='SYMBOL', show=False, return_fig=True,
            save=f"-curated_celltype_markers_quantile-{sample_name}-{cell_type}.png")
        # move output file to correct location
        if not os.path.exists(snakemake.output["dir"]):
            os.mkdir(snakemake.output["dir"])
        os.rename(
            f"show-curated_celltype_markers_quantile-{sample_name}-{cell_type}.png",
            snakemake.output["dir"] + f"/{sample_name}-{cell_type}.png")
