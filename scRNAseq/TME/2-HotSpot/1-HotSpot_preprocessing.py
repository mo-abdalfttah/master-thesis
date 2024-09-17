# Load Libraries

import scanpy as sc
import scipy
import mplscience
import hotspot
import pickle
import matplotlib.pyplot as scplt
import anndata as ann
import os

# Set Working Directory

os.chdir("/home/groups/singlecell/mabdalfttah/projects/Master's_Thesis/")
path_dir = "scRNAseq/TME/2-HotSpot/"

# Read Dataset
adata_norm = sc.read_h5ad("scRNAseq/TME/1-QC/R_objects_2023-11-24/MSS_Pellka.h5ad")
adata_norm
# Check point since we converted this adata object from seurat
adata_norm.obsm['PCA'] 
adata_norm.obsm['X_pca'] = adata_norm.obsm['PCA'].values 
sc.pl.pca(adata_norm, color='CST3')
adata_norm.obsm['X_umap'] = adata_norm.obsm['UMAP'].values 
sc.pl.umap(adata_norm, color=['SampleID'])
sc.pl.umap(adata_norm, color=['MMRStatus'])

adata_norm.layers['counts'] = adata_norm.X.copy()
adata_norm.layers["counts_csc"] = scipy.sparse.csc_matrix(adata_norm.layers["counts"])

### Chekc we use the counts not the normalized
adata_norm.to_df()
adata_norm.to_df(layer="logcounts")
# Create patient files
set(adata_norm.obs['Patient'])

# Create Directories for each Patient
import os
for patient in set(adata_norm.obs['Patient']):
    path = path_dir + patient + "/"
    try:
        os.mkdir(path)
    except OSError as error:
        print(error) 


# Split the object by Paients and save it
for patient in set(adata_norm.obs['Patient']):
    currdata = adata_norm[adata_norm.obs['Patient'] == patient]
    currdata.obsm['X_pca'].shape
    sc.pp.filter_genes(currdata, min_counts=1)
    d = list([x for x in currdata.var_names if list(currdata.var_names).count(x) > 1])
    keep_genes = set(currdata.var_names)-set(d)
    currdata = currdata[:,currdata.var_names.isin(list(keep_genes))]
    currdata
    filename = path_dir + patient + "/" + patient + ".h5ad"
    currdata.write(filename, compression='gzip')



# Run the Preprocessing and save the PCA

for patient in set(adata_norm.obs['Patient']):
    filename = path_dir + patient + "/" + patient + ".h5ad"
    currdata = ann.read_h5ad(filename)
    currdata
    sc.pp.normalize_total(currdata)
    sc.pp.log1p(currdata)
    sc.pp.scale(currdata)
    sc.tl.pca(currdata, n_comps=10)
    filename = "PCA" + patient + "_MSS.pdf"
    with mplscience.style_context():
        sc.pl.pca_variance_ratio(currdata, save= filename)

# Create and Save Hotspot object
for patient in set(adata_norm.obs['Patient']):
    adata = ann.read_h5ad((path_dir+patient+"/"+patient+".h5ad"))
    adata
    hs = hotspot.Hotspot(
        adata,
        layer_key='counts_csc',
        model='danb',
        latent_obsm_key="X_pca",
        umi_counts_obs_key="nCount_RNA")
    hs
    filename = path_dir + patient + "/" + patient + ".p"
    with open(filename, 'wb') as file:
        pickle.dump(hs, file)


# Compute Knn Graph
for patient in set(adata_norm.obs['Patient']):
    filename = path_dir + patient + "/" + patient + ".p"
    with open(path_dir+patient+"/"+patient+".p", 'rb') as f:
        hs = pickle.load(f)
    hs
    hs.create_knn_graph(weighted_graph=False, n_neighbors= 30)
    hs.neighbors
    with open(filename, 'wb') as file:
        pickle.dump(hs, file)

# Determining informative genes
for patient in set(adata_norm.obs['Patient']):
    filename = path_dir + patient + "/" + patient + ".p"
    with open(path_dir+ patient + "/" + patient + ".p", 'rb') as f:
        hs = pickle.load(f)
    hs
    hs_results = hs.compute_autocorrelations(jobs=40)
    hs_results.to_csv(path_dir + patient + "/" + "Autocorrelations_" + patient + "_30knn.csv")
    with open(filename, 'wb') as file:
        pickle.dump(hs, file)
