import pandas as pd
import numpy as np
import SpatialDE
import NaiveDE
import os

# Weird spatialDE ravel fix.
if not hasattr(pd.Series, "ravel"):
    pd.Series.ravel = lambda self: self.to_numpy().ravel()

# Define the output directory and ensure it exists globally
output_dir = "SpatialDE_Results"
os.makedirs(output_dir, exist_ok=True)

datasets = {
    "Rep11_MOB": {
        "type": "tsv_counts",
        "file": "data/Rep11_MOB_count_matrix-1.tsv",
    },
    "Layer_2_BC": {
        "type": "tsv_counts", 
        "file": "data/Layer2_BC_count_matrix.tsv", 
    },
    "seqfish": {
        "type": "xlsx", 
        "file": "data/seqfish.xlsx"
    },
}

for dataset_name, info in datasets.items():
    print(f"\n{'='*50}\nProcessing dataset: {dataset_name}\n{'='*50}")

    file_unscaled = f"{output_dir}/spatialDE_{dataset_name}_unscaled_results.csv"
    file_scaled = f"{output_dir}/spatialDE_{dataset_name}_scaled_results.csv"
    
    if os.path.exists(file_unscaled) and os.path.exists(file_scaled):
        print(f"Skipping {dataset_name}: Both scaled and unscaled results already exist.")
        continue

    # ==========================================
    # 1. LOAD DATA 
    # ==========================================
    if info["type"] == "tsv_counts":
        if not os.path.exists(info["file"]):
            print(f"Skipping {dataset_name}: File not found.")
            continue

        # Load raw expression (Index 0 is the Spot IDs like 16.92x9.015)
        df = pd.read_csv(info["file"], sep="\t", index_col=0)
        
        # NO TRANSPOSE NEEDED: The raw file is already Spots (rows) x Genes (columns)
        counts = df 

        # Extract coordinates by splitting the spot names (e.g. "16.92x9.015")
        coords = counts.index.to_series().str.split('x', expand=True)
        sample_info = pd.DataFrame({
            'x': coords[0].astype(float),
            'y': coords[1].astype(float)
        }, index=counts.index)

    elif info["type"] == "xlsx":
        if not os.path.exists(info["file"]):
            print(f"Skipping {dataset_name}: File not found.")
            continue

        counts_df = pd.read_excel(info["file"], sheet_name="Hippocampus Counts", header=None)
        counts_df.rename(columns={0: "gene"}, inplace=True)
        counts_df.set_index("gene", inplace=True)
        counts = counts_df.T

        coord_df = pd.read_excel(info["file"], sheet_name="Centroids", header=None)
        sample_info = coord_df.iloc[:, :2].copy()
        sample_info.columns = ["x", "y"]
        sample_info.index = counts.index

    else:
        print(f"Unknown file type for {dataset_name}. Skipping.")
        continue

    # ==========================================
    # 2. FILTER AND NORMALIZE
    # ==========================================
    counts = counts.loc[:, counts.sum(axis=0) >= 3]
    
    sample_info["total_counts"] = counts.sum(axis=1)
    
    valid_spots = sample_info["total_counts"] > 10
    sample_info = sample_info[valid_spots]
    counts = counts.loc[valid_spots]
    
    print(f"Filtered matrix shape: {counts.shape[0]} cells/spots, {counts.shape[1]} genes.")

    sample_info["x"] = sample_info["x"].astype(float)
    sample_info["y"] = sample_info["y"].astype(float)
    sample_info = sample_info.loc[counts.index]

    # Takes the raw counts SPARK used and stabilizes/regresses them for SpatialDE.
    print("Normalizing and regressing out sequencing depth...")
    norm_counts = NaiveDE.stabilize(counts)
    res = NaiveDE.regress_out(sample_info, norm_counts.T, "np.log(total_counts)").T

    # ==========================================
    # 3. RUN SPATIALDE (SCALED VS UNSCALED)
    # ==========================================
    for is_scaled in [False, True]:
        scale_label = "scaled" if is_scaled else "unscaled"
        output_filename = f"{output_dir}/spatialDE_{dataset_name}_{scale_label}_results.csv"

        if os.path.exists(output_filename):
            print(f"\n-> Skipping SpatialDE ({scale_label}): File already exists.")
            continue

        print(f"\n--- Running SpatialDE for {dataset_name} ({scale_label}) ---")
        current_sample_info = sample_info.copy()

        if is_scaled:
            print("Scaling coordinates to [0, 1] range...")
            current_sample_info["x"] = (current_sample_info["x"] - current_sample_info["x"].min()) / (
                current_sample_info["x"].max() - current_sample_info["x"].min()
            )
            current_sample_info["y"] = (current_sample_info["y"] - current_sample_info["y"].min()) / (
                current_sample_info["y"].max() - current_sample_info["y"].min()
            )
        else:
            print("Using raw coordinates...")

        X = current_sample_info[["x", "y"]].values

        results = SpatialDE.run(X, res)

        results_sorted = results.sort_values("qval")
        print(f"\nTop 5 Spatially Variable Genes ({scale_label}):")
        print(results_sorted[["g", "l", "pval", "qval"]].head(5))

        results_sorted.to_csv(output_filename, index=False)
        print(f"Saved results to {output_filename}\n")