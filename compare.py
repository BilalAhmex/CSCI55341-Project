import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2
import os
import random
from scipy.spatial.distance import pdist, squareform

# ── Custom Benchmarking Functions ──────────────────────────────────────────

def rank_biased_overlap(list1, list2, p=0.9):
    s1, s2 = set(), set()
    score = 0.0
    max_depth = max(len(list1), len(list2))
    
    for d in range(1, max_depth + 1):
        if d - 1 < len(list1): s1.add(list1[d - 1])
        if d - 1 < len(list2): s2.add(list2[d - 1])
        agreement = len(s1.intersection(s2)) / d
        score += (p**(d - 1)) * agreement
        
    return (1 - p) * score

def calculate_morans_i(expr, coords):
    n = len(expr)
    dists = squareform(pdist(coords))
    np.fill_diagonal(dists, np.inf)
    w = 1.0 / dists
    np.fill_diagonal(w, 0)
    
    row_sums = w.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1 
    w = w / row_sums
    
    W_sum = w.sum()
    x_bar = np.mean(expr)
    diff = expr - x_bar
    
    numerator = np.sum(w * np.outer(diff, diff))
    denominator = np.sum(diff**2)
    
    if denominator == 0:
        return 0.0
        
    return (n / W_sum) * (numerator / denominator)

# ── Load results ──────────────────────────────────────────────────────────────
spark_file = "SPARK_Results/SPARK_Rep11_MOB_scaled_results.csv"
spatialde_file = "SpatialDE_Results/spatialDE_Rep11_MOB_scaled_results.csv"

spark_df      = pd.read_csv(spark_file)
spatialde_df  = pd.read_csv(spatialde_file)

spatialde_df = spatialde_df.rename(columns={"g": "gene", "pval": "pval_spatialde", "qval": "qval_spatialde"})

if "combined_pvalue" in spark_df.columns:
    spark_df = spark_df.rename(columns={"combined_pvalue": "pval_spark", "adjusted_pvalue": "qval_spark"})

# ── Significant gene sets ─────────────────────────────────────────────────────
spark_sig     = set(spark_df[spark_df["qval_spark"]     < 0.05]["gene"])
spatialde_sig = set(spatialde_df[spatialde_df["qval_spatialde"] < 0.05]["gene"])

both           = spark_sig & spatialde_sig
only_spark     = spark_sig - spatialde_sig
only_spatialde = spatialde_sig - spark_sig

print(f"Significant in SPARK only:     {len(only_spark)}")
print(f"Significant in SpatialDE only: {len(only_spatialde)}")
print(f"Significant in both:           {len(both)}")

# ── Merge side by side ────────────────────────────────────────────────────────
merged = pd.merge(
    spark_df[["gene", "pval_spark", "qval_spark"]],
    spatialde_df[["gene", "pval_spatialde", "qval_spatialde", "FSV", "l"]],
    on  = "gene",
    how = "outer"
)

merged["sig_spark"]     = merged["qval_spark"]     < 0.05
merged["sig_spatialde"] = merged["qval_spatialde"] < 0.05
merged["agreement"]     = merged["sig_spark"] == merged["sig_spatialde"]

# ── Benchmarking 1: Rank-Biased Overlap ──────────────────────────────────────
spark_ranked = spark_df.sort_values("qval_spark")["gene"].tolist()
spatialde_ranked = spatialde_df.sort_values("qval_spatialde")["gene"].tolist()

rbo_score = rank_biased_overlap(spark_ranked[:500], spatialde_ranked[:500], p=0.95)
print(f"\nRank-Biased Overlap (Top 500): {rbo_score:.4f}")

common = merged.dropna(subset=["qval_spark", "qval_spatialde"])
corr = common["qval_spark"].rank().corr(common["qval_spatialde"].rank(), method="spearman")
print(f"Spearman correlation of rankings: {corr:.4f}")

# ── Benchmarking 3: Moran's I Data Collection ─────────────────────────────────
print("\n--- Calculating Moran's I for Boxplot Comparison ---")
expr_file = "data/Rep11_MOB_trans.tsv"
coord_file = "data/Rep11_MOB_trans.idx"

mi_spark, mi_spatialde, mi_both = [], [], []

if os.path.exists(expr_file) and os.path.exists(coord_file):
    counts_df = pd.read_csv(expr_file, sep='\t').set_index('gene')
    if 'ensemblid' in counts_df.columns:
        counts_df = counts_df.drop('ensemblid', axis=1)
    
    coord_df = pd.read_csv(coord_file, sep='\t', index_col=0)
    coords = coord_df.T[['x', 'y']].astype(float).values
    
    # Sample up to 30 genes from each group to keep the script fast
    random.seed(42)
    n_samp = 30
    samp_spark = random.sample(list(only_spark), min(len(only_spark), n_samp))
    samp_spatialde = random.sample(list(only_spatialde), min(len(only_spatialde), n_samp))
    samp_both = random.sample(list(both), min(len(both), n_samp))
    
    for g in samp_spark:
        if g in counts_df.index: mi_spark.append(calculate_morans_i(counts_df.loc[g].values, coords))
    for g in samp_spatialde:
        if g in counts_df.index: mi_spatialde.append(calculate_morans_i(counts_df.loc[g].values, coords))
    for g in samp_both:
        if g in counts_df.index: mi_both.append(calculate_morans_i(counts_df.loc[g].values, coords))
else:
    print("Raw data not found. Skipping Moran's I step.")

# ── Color coding & P-value Safety ─────────────────────────────────────────────
def assign_color(row):
    if row["sig_spark"] and row["sig_spatialde"]: return "darkgreen"
    elif row["sig_spark"]: return "blue"
    elif row["sig_spatialde"]: return "red"
    else: return "lightgrey"

common = common.copy()
common["color"] = common.apply(assign_color, axis=1)

# Prevent log10(0) by clipping to the smallest possible python float
eps = 1e-300
common["safe_qval_spark"] = np.clip(common["qval_spark"], eps, 1.0)
common["safe_qval_spatialde"] = np.clip(common["qval_spatialde"], eps, 1.0)

# ── CREATE OUTPUT DIRECTORY ───────────────────────────────────────────────────
comp_dir = "Comparison"
os.makedirs(comp_dir, exist_ok=True)

# ── INDIVIDUAL PLOTTING ───────────────────────────────────────────────────────

# Plot 1: P-Value Scatter Plot
plt.figure(figsize=(8, 6))
plt.scatter(
    -np.log10(common["safe_qval_spatialde"]),
    -np.log10(common["safe_qval_spark"]),
    c=common["color"], alpha=0.7, s=30, edgecolors="none"
)
plt.axhline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
plt.axvline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
plt.xlabel("-log10 SpatialDE qval")
plt.ylabel("-log10 SPARK adjusted p")
plt.title("SPARK vs SpatialDE SVG Agreement")

legend_patches = [
    mpatches.Patch(color="darkgreen", label=f"Both ({len(both)})"),
    mpatches.Patch(color="blue",      label=f"SPARK only ({len(only_spark)})"),
    mpatches.Patch(color="red",       label=f"SpatialDE only ({len(only_spatialde)})"),
]
plt.legend(handles=legend_patches, fontsize=9)

top_both = common[common["color"] == "darkgreen"].nsmallest(5, "qval_spark")
for _, row in top_both.iterrows():
    plt.annotate(row["gene"], xy=(-np.log10(row["safe_qval_spatialde"]), -np.log10(row["safe_qval_spark"])),
                 fontsize=8, color="black", xytext=(5, 5), textcoords="offset points")
plt.tight_layout()
scatter_path = os.path.join(comp_dir, "SPARK_vs_SpatialDE_Scatter.png")
plt.savefig(scatter_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved Scatter Plot: {scatter_path}")

# Plot 2: Venn Diagram
plt.figure(figsize=(6, 6))
venn2(subsets=[spark_sig, spatialde_sig], set_labels=("SPARK", "SpatialDE"))
plt.title(f"SVG Overlap (qval < 0.05)\nRBO Score: {rbo_score:.4f}")
plt.tight_layout()
venn_path = os.path.join(comp_dir, "SPARK_vs_SpatialDE_Venn.png")
plt.savefig(venn_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved Venn Diagram: {venn_path}")

# Plot 3: Moran's I Boxplot
if mi_spark or mi_spatialde or mi_both:
    plt.figure(figsize=(7, 6))
    data = [mi_spark, mi_spatialde, mi_both]
    labels = [f"SPARK Only\n(n={len(mi_spark)})", f"SpatialDE Only\n(n={len(mi_spatialde)})", f"Both\n(n={len(mi_both)})"]
    
    bplot = plt.boxplot(data, patch_artist=True, labels=labels)
    colors = ['blue', 'red', 'darkgreen']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        
    plt.ylabel("Moran's I Score")
    plt.title("Spatial Autocorrelation by Group")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    moran_path = os.path.join(comp_dir, "SPARK_vs_SpatialDE_MoransI.png")
    plt.savefig(moran_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved Moran's I Boxplot: {moran_path}")

# ── Save merged comparison CSV ────────────────────────────────────────────────
merged_sorted = merged.sort_values("qval_spark")
csv_path = os.path.join(comp_dir, "SPARK_vs_SpatialDE_merged.csv")
merged_sorted.to_csv(csv_path, index=False)
print(f"\nMerged results saved to {csv_path}")