import scanpy as sc
import scipy
import mplscience
import hotspot
import pickle
import matplotlib.pyplot as plt
import anndata as ann
import os

os.chdir("/home/groups/singlecell/mabdalfttah/projects/Master's_Thesis/")
path_dir = "scRNAseq/TME/2-HotSpot/"

adata_norm = sc.read_h5ad("scRNAseq/TME/1-QC/R_objects_2023-11-24/MSS_Pellka.h5ad")

patients = ['C171']

# for patient in set(adata_norm.obs['Patient'])


for patient in set(patients):
    filename = path_dir + patient + "/" + patient + ".p"
    with open(path_dir+patient+"/"+patient+".p", 'rb') as f:
        hs = pickle.load(f)
    hs
    adata = ann.read_h5ad((path_dir+patient+"/"+patient+".h5ad"))
    adata
    for param in [200,150,100,50]:
        modules = hs.create_modules(min_gene_threshold=param, core_only=False, fdr_threshold=0.05)
        modules.value_counts()
        modules.to_csv((path_dir+patient+"/Modules_"+patient+"_CORE_30knn_"+str(param)+"genes.csv"))
        hs.plot_local_correlations(vmin=-12, vmax=12)
        plt.savefig((path_dir + "Hotspot_Heatmap_"+patient+"_CORE_30knn_"+str(param)+"genes.pdf"))
        module_scores = hs.calculate_module_scores()
        module_scores.to_csv((path_dir+patient+"/ModuleScores_"+patient+"_CORE_30knn_"+str(param)+"genes.csv"))
        module_scores = hs.module_scores
        module_cols = []
        for c in module_scores.columns:
            key = f"Module {c}"
            adata.obs[key] = module_scores[c]
            module_cols.append(key)
            with mplscience.style_context():
                sc.pl.umap(adata, color=module_cols, frameon=False, vmin=-12, vmax=12, save=("UMAP_ModuleScores_"+patient+"_CORE_30knn_"+str(param)+"genes.pdf"))
    adata.write(path_dir+patient+"/"+patient+".h5ad", compression='gzip')
    with open(filename, 'wb') as file:
        pickle.dump(hs, file)


import os
