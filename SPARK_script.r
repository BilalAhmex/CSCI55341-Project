library(SPARK)
library(Matrix)

# ── File paths ────────────────────────────────────────────────────────────────
expression_file <- "data/Rep11_MOB_trans.tsv"
coordinate_file <- "data/Rep11_MOB_trans.idx"

# ── Load expression data ──────────────────────────────────────────────────────
df <- read.delim(expression_file, sep = "\t", row.names = "gene", check.names = FALSE)

# Drop ensemblid column if present
if ("ensemblid" %in% colnames(df)) {
  df <- df[, colnames(df) != "ensemblid"]
}

# Transpose: spots as rows, genes as columns (same as your counts = df.T)
counts <- t(df)  # now: spots x genes

# ── Filter lowly expressed genes (sum >= 3) ───────────────────────────────────
gene_sums <- colSums(counts)
counts <- counts[, gene_sums >= 3]
cat(sprintf("Filtered matrix: %d spots, %d genes.\n", nrow(counts), ncol(counts)))

# ── Load and process coordinates ──────────────────────────────────────────────
coord_df <- read.delim(coordinate_file, sep = "\t", row.names = 1, check.names = FALSE)

# Transpose so spots are rows, x and y are columns (same as sample_info = coord_df.T)
sample_info <- as.data.frame(t(coord_df))
sample_info$x <- as.numeric(sample_info$x)
sample_info$y <- as.numeric(sample_info$y)

# Scale x and y to [0, 1] range
sample_info$x <- (sample_info$x - min(sample_info$x)) / (max(sample_info$x) - min(sample_info$x))
sample_info$y <- (sample_info$y - min(sample_info$y)) / (max(sample_info$y) - min(sample_info$y))

# Align coordinates with expression matrix rows
sample_info <- sample_info[rownames(counts), ]

# ── Prepare count matrix for SPARK ───────────────────────────────────────────
# SPARK expects genes x spots (opposite of your counts matrix), so transpose back
counts_spark <- t(counts)  # now: genes x spots

# ── Create SPARK object ───────────────────────────────────────────────────────
spark <- CreateSPARKObject(
  counts      = counts_spark,
  location    = sample_info[, c("x", "y")],
  percentage  = 0.1,           # min fraction of spots expressing a gene
  min_total_counts = 3         # matches your sum >= 3 filter
)

spark@lib_size <- apply(spark@counts, 2, sum)

cat(sprintf("SPARK object: %d genes, %d spots.\n",
            nrow(spark@counts), ncol(spark@counts)))

# ── Run SPARK ─────────────────────────────────────────────────────────────────
cat("Fitting variance components...\n")
spark <- spark.vc(
  spark,
  covariates       = NULL,
  lib_size         = spark@lib_size,
  num_core         = 4,       # adjust to your machine
  verbose          = FALSE
)

cat("Running hypothesis tests...\n")
spark <- spark.test(
  spark,
  check_positive = TRUE,
  verbose        = FALSE
)

# ── View and save results ─────────────────────────────────────────────────────
results <- spark@res_mtest
results$gene <- rownames(results)

# Sort by adjusted p-value (equivalent to qval in SpatialDE)
results_sorted <- results[order(results$adjusted_pvalue), ]

cat("\nTop 5 Spatially Variable Genes:\n")
print(head(results_sorted[, c("gene", "combined_pvalue", "adjusted_pvalue")], 5))

write.csv(results_sorted, "SPARK_Rep11_MOB_results.csv", row.names = FALSE)
cat("Results saved to SPARK_Rep11_MOB_results.csv\n")