##
# imports
##

import pandas as pd

##
# Import fgsea results
##

fgsea = pd.read_csv(snakemake.input['fgsea'], sep="\t")
# fgsea.head()

##
# Filter fgsea results
##

search_patterns = [
    "vascular",
    "endothelial",
    "fibroblast",
    "skeletal",
    "muscle",
    "immune",
    "macrophage",
    # "stromal"
]

fgsea['pattern_hit'] = [any(pattern in pathway.lower() for pattern in search_patterns) for pathway in fgsea['pathway']]
# fgsea.head()

top_hits = fgsea[fgsea['pattern_hit'] == True].drop_duplicates(['group'])
# top_hits.head()

##
# Produce relabelled spatial plot
##

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
genes = pd.read_csv(snakemake.input['genes'], sep="\t")

##
# import information about samples
##

samples = (
    pd.read_csv("config/samples.tsv", sep="\t", dtype={"sample_name": str})
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
sc.pp.neighbors(slide, n_neighbors=10, n_pcs=40)
sc.tl.umap(slide)
sc.tl.leiden(slide, resolution=samples['resolution'][sample_name], key_added="clusters")
sc.tl.rank_genes_groups(slide, 'clusters', method='t-test')

def label_cluster(clusters, mapping_table):
    mapping_dict = {}
    for i in range(len(mapping_table)):
        mapping_dict[str(list(mapping_table['group'])[i])] = list(mapping_table['pathway'])[i]
    out = []
    for cluster in clusters:
        if cluster in mapping_dict.keys():
            out.append(mapping_dict[cluster])
        else:
            out.append(cluster)
    # out = pd.Series(out)
    return out

slide.obs.insert(slide.obs.shape[1], 'labels', label_cluster(slide.obs['clusters'], top_hits))
# slide.obs.head()

##
# Produce plots
##

with mpl.rc_context({'figure.figsize': [12,6],
                     'axes.facecolor': 'white',
                     "savefig.bbox": 'tight'}):
    fig = sc.pl.spatial(slide, img_key="hires", color="labels", size=1.5)
    plt.savefig(snakemake.output['png'])
