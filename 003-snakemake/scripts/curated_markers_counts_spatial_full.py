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

def read_and_qc(sample_name, s_col='sample'):
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
    
    s_keys = list(adata.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[sample_name in k for k in s_keys]][0]
    
    adata.uns['spatial'] = {s_spatial: adata.uns['spatial'][s_spatial]}
    
    return adata

##
# Load curated markers
##

curated_markers = pd.read_csv(snakemake.input["tsv"], sep="\t")
curated_markers.head()

##
# Process each sample
##

for sample_name in sample_names:
    slide = read_and_qc(sample_name)
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
            vmin=0, vmax=None,
            gene_symbols='SYMBOL', show=False, return_fig=True,
            save=f"-curated_celltype_markers_full-{sample_name}-{cell_type}.png")
        # move output file to correct location
        if not os.path.exists(snakemake.output["dir"]):
            os.mkdir(snakemake.output["dir"])
        os.rename(
            f"show-curated_celltype_markers_full-{sample_name}-{cell_type}.png",
            snakemake.output["dir"] + f"/{sample_name}-{cell_type}.png")
