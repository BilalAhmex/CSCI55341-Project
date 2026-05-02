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
        "type": "tsv",
        "expr_file": "data/Rep11_MOB_trans.tsv",
        "coord_file": "data/Rep11_MOB_trans.idx",
    },
    "Layer_2_BC": {
        "type": "tsv",
        "expr_file": "data/Layer2_BC_trans.tsv",
        "coord_file": "data/Layer2_BC_trans.idx",
    },
    "seqfish": {
        "type": "xlsx", 
        "file": "data/seqfish.xlsx"
    },
}

for dataset_name, info in datasets.items():
    print(f"\n{'='*50}\nProcessing dataset: {dataset_name}\n{'='*50}")

    # Check if both scaled and unscaled results already exist to save time
    file_unscaled = f"{output_dir}/spatialDE_{dataset_name}_unscaled_results.csv"
    file_scaled = f"{output_dir}/spatialDE_{dataset_name}_scaled_results.csv"
    
    if os.path.exists(file_unscaled) and os.path.exists(file_scaled):
        print(f"Skipping {dataset_name}: Both scaled and unscaled results already exist.")
        continue

    # ==========================================
    # 1. LOAD DATA
    # ==========================================
    if info["type"] == "tsv":
        if not os.path.exists(info["expr_file"]) or not os.path.exists(
            info["coord_file"]
        ):
            print(f"Skipping {dataset_name}: Files not found.")
            continue

        # Load expression
        df = pd.read_csv(info["expr_file"], sep="\t")
        df = df.set_index("gene")
        if "ensemblid" in df.columns:
            df = df.drop("ensemblid", axis=1)
        counts = df.T

        # Load coordinates
        coord_df = pd.read_csv(info["coord_file"], sep="\t", index_col=0)
        sample_info = coord_df.T

    elif info["type"] == "xlsx":
        if not os.path.exists(info["file"]):
            print(f"Skipping {dataset_name}: File not found.")
            continue

        # Load expression
        counts_df = pd.read_excel(
            info["file"], sheet_name="Hippocampus Counts", header=None
        )
        counts_df.rename(columns={0: "gene"}, inplace=True)
        counts_df.set_index("gene", inplace=True)
        counts = counts_df.T

        # Load coordinates
        coord_df = pd.read_excel(info["file"], sheet_name="Centroids", header=None)

        # Grab only the first two columns and name them x and y
        sample_info = coord_df.iloc[:, :2].copy()
        sample_info.columns = ["x", "y"]

        # Match the integer indices from the counts matrix to the coordinates
        sample_info.index = counts.index

    else:
        print(f"Unknown file type for {dataset_name}. Skipping.")
        continue

    # ==========================================
    # 2. FILTER AND NORMALIZE
    # ==========================================
    # Filter out lowly expressed genes for speed
    counts = counts.loc[:, counts.sum(axis=0) >= 3]
    
    # Calculate total counts per cell (sequencing depth)
    sample_info["total_counts"] = counts.sum(axis=1)
    
    # Filter out empty or extremely poor-quality spots
    valid_spots = sample_info["total_counts"] > 10
    sample_info = sample_info[valid_spots]
    counts = counts.loc[valid_spots]
    
    print(
        f"Filtered matrix shape: {counts.shape[0]} cells/spots, {counts.shape[1]} genes."
    )

    # Ensure the data types are floats
    sample_info["x"] = sample_info["x"].astype(float)
    sample_info["y"] = sample_info["y"].astype(float)

    # Ensure the coordinates perfectly align with the expression matrix rows
    sample_info = sample_info.loc[counts.index]

    # Normalize and regress out sequencing depth using NaiveDE
    # (We only need to do this math once per dataset!)
    print("Normalizing and regressing out sequencing depth...")
    norm_counts = NaiveDE.stabilize(counts)
    res = NaiveDE.regress_out(sample_info, norm_counts.T, "np.log(total_counts)").T


    # ==========================================
    # 3. RUN SPATIALDE (SCALED VS UNSCALED)
    # ==========================================
    for is_scaled in [False, True]:
        scale_label = "scaled" if is_scaled else "unscaled"
        output_filename = f"{output_dir}/spatialDE_{dataset_name}_{scale_label}_results.csv"

        # Check if this specific configuration already exists
        if os.path.exists(output_filename):
            print(f"\n-> Skipping SpatialDE ({scale_label}): File already exists.")
            continue

        print(f"\n--- Running SpatialDE for {dataset_name} ({scale_label}) ---")
        
        # Create a copy of the coordinates so we don't permanently alter them
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

        # Extract spatial coordinates for this run
        X = current_sample_info[["x", "y"]].values

        # Run SpatialDE on the regressed data
        results = SpatialDE.run(X, res)

        # View results
        results_sorted = results.sort_values("qval")
        print(f"\nTop 5 Spatially Variable Genes ({scale_label}):")
        print(results_sorted[["g", "l", "pval", "qval"]].head(5))

        # Save the results to folder
        results_sorted.to_csv(output_filename, index=False)
        print(f"Saved results to {output_filename}\n")