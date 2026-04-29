import pandas as pd
import numpy as np
import SpatialDE

# Weird spatialDE ravel fix.
if not hasattr(pd.Series, 'ravel'):
    pd.Series.ravel = lambda self: self.to_numpy().ravel()


# File paths
expression_file = 'data/Rep11_MOB_trans.tsv' 
coordinate_file = 'data/Rep11_MOB_trans.idx'

# Load and clean
df = pd.read_csv(expression_file, sep='\t')

# Set gene as the row names and remove the ensemblid column
df = df.set_index('gene')
if 'ensemblid' in df.columns:
    df = df.drop('ensemblid', axis=1)

# Transpose so spots are rows and genes are columns (spatialDE requirment)
counts = df.T 

# Filter out lowly expressed genes for speed
counts = counts.loc[:, counts.sum(axis=0) >= 3]
print(f"Filtered matrix shape: {counts.shape[0]} spots, {counts.shape[1]} genes.")


# Read the .idx file, setting the first column ('coordinate', 'x', 'y') as the index
coord_df = pd.read_csv(coordinate_file, sep='\t', index_col=0)

# Transpose it so Spots are Rows, and 'x' and 'y' become the columns
sample_info = coord_df.T

# Ensure the data types are floats and not strings
sample_info['x'] = sample_info['x'].astype(float)
sample_info['y'] = sample_info['y'].astype(float)

# Scale X and Y coordinates to the [0, 1] range
sample_info['x'] = (sample_info['x'] - sample_info['x'].min()) / (sample_info['x'].max() - sample_info['x'].min())
sample_info['y'] = (sample_info['y'] - sample_info['y'].min()) / (sample_info['y'].max() - sample_info['y'].min())

# Ensure the coordinates perfectly align with the expression matrix rows
sample_info = sample_info.loc[counts.index]

# Run SpatialDE
X = sample_info[['x', 'y']].values

print("Running SpatialDE: ")
# Since SpatialDB's values are already decimals (normalized), we feed counts directly into SpatialDE without running NaiveDE for normalization
results = SpatialDE.run(X, counts)

# View and save results
results_sorted = results.sort_values('qval')

print("\nTop 5 Spatially Variable Genes:")
print(results_sorted[['g', 'l', 'pval', 'qval']].head(5))

# Save the results
results_sorted.to_csv('spatialDE_Rep11_MOB_results.csv', index=False)