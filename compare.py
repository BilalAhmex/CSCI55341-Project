import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2

# ── Load results ──────────────────────────────────────────────────────────────
spark_df      = pd.read_csv("SPARK_Rep11_MOB_results.csv")
spatialde_df  = pd.read_csv("spatialDE_Rep1_MOB_results.csv")

# Rename SpatialDE columns for clarity
spatialde_df = spatialde_df.rename(columns={
    "g"    : "gene",
    "pval" : "pval_spatialde",
    "qval" : "qval_spatialde"
})

# Rename SPARK columns for clarity
spark_df = spark_df.rename(columns={
    "combined_pvalue"  : "pval_spark",
    "adjusted_pvalue"  : "qval_spark"
})

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

print("\nTop 10 by SPARK adjusted p-value:")
print(merged.sort_values("qval_spark").head(10)[
    ["gene", "qval_spark", "qval_spatialde", "FSV", "l"]
])

# ── Spearman correlation of rankings ─────────────────────────────────────────
common = merged.dropna(subset=["qval_spark", "qval_spatialde"])
corr = common["qval_spark"].rank().corr(common["qval_spatialde"].rank(), method="spearman")
print(f"\nSpearman correlation of p-value rankings: {corr:.4f}")

# ── Color coding for plots ────────────────────────────────────────────────────
def assign_color(row):
    if row["sig_spark"] and row["sig_spatialde"]:
        return "darkgreen"   # both
    elif row["sig_spark"]:
        return "blue"        # SPARK only
    elif row["sig_spatialde"]:
        return "red"         # SpatialDE only
    else:
        return "lightgrey"   # neither

common = common.copy()
common["color"] = common.apply(assign_color, axis=1)

# ── Plot 1: p-value scatter ───────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
ax.scatter(
    -np.log10(common["qval_spatialde"]),
    -np.log10(common["qval_spark"]),
    c     = common["color"],
    alpha = 0.7,
    s     = 30,
    edgecolors = "none"
)
ax.axhline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
ax.axvline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
ax.set_xlabel("-log10 SpatialDE qval")
ax.set_ylabel("-log10 SPARK adjusted p")
ax.set_title("SPARK vs SpatialDE SVG Agreement")

legend_patches = [
    mpatches.Patch(color="darkgreen", label=f"Both ({len(both)})"),
    mpatches.Patch(color="blue",      label=f"SPARK only ({len(only_spark)})"),
    mpatches.Patch(color="red",       label=f"SpatialDE only ({len(only_spatialde)})"),
    mpatches.Patch(color="lightgrey", label="Neither"),
]
ax.legend(handles=legend_patches, fontsize=8)

# Annotate top 5 genes agreed upon by both
top_both = common[common["color"] == "darkgreen"].nsmallest(5, "qval_spark")
for _, row in top_both.iterrows():
    ax.annotate(
        row["gene"],
        xy        = (-np.log10(row["qval_spatialde"]), -np.log10(row["qval_spark"])),
        fontsize  = 7,
        color     = "black",
        xytext    = (5, 5),
        textcoords= "offset points"
    )

# ── Plot 2: Venn diagram ──────────────────────────────────────────────────────
ax2 = axes[1]
venn2(
    subsets = [spark_sig, spatialde_sig],
    set_labels = ("SPARK", "SpatialDE"),
    ax = ax2
)
ax2.set_title("SVG Overlap (qval < 0.05)")

plt.tight_layout()
plt.savefig("SPARK_vs_SpatialDE_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("Plot saved to SPARK_vs_SpatialDE_comparison.png")

# ── Save merged comparison CSV ────────────────────────────────────────────────
merged_sorted = merged.sort_values("qval_spark")
merged_sorted.to_csv("SPARK_vs_SpatialDE_merged.csv", index=False)
print("Merged results saved to SPARK_vs_SpatialDE_merged.csv")