# run_analysis.py
import sys
import scanpy as sc
import mplscience
import hotspot
import pickle
import os

os.chdir("/home/groups/singlecell/mabdalfttah/projects/Master's_Thesis/")
path_dir = "scRNAseq/TME/2-HotSpot/"

def process_patient(patient):
    filename = path_dir + patient + "/" + patient + ".p"
    with open(path_dir+patient+"/"+patient+".p", 'rb') as f:
        hs = pickle.load(f)

    hs_results = hs.results
    hs_genes = hs_results.loc[hs_results.FDR < 0.05].index
    local_correlations = hs.compute_local_correlations(hs_genes, jobs=40)
    local_correlations.to_csv(path_dir+patient+"/"+"Localcorrelations_"+patient+"_30knn.csv")

    with open(filename, 'wb') as file:
        pickle.dump(hs, file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_analysis.py <patient_name>")
        sys.exit(1)

    patient_name = sys.argv[1]
    process_patient(patient_name)
