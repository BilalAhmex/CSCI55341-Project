library(SPARK)
library(Matrix)

# ── Output directory ──────────────────────────────────────────────────────────
output_dir <- "SPARK_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Dataset definitions ───────────────────────────────────────────────────────
datasets <- list(
  Rep11_MOB = list(
    type       = "tsv",
    expr_file  = "data/Rep11_MOB_trans.tsv",
    coord_file = "data/Rep11_MOB_trans.idx"
  ),
  Layer_2_BC = list(
    type       = "tsv",
    expr_file  = "data/Layer2_BC_trans.tsv",
    coord_file = "data/Layer2_BC_trans.idx"
  ),
  seqfish = list(
    type = "xlsx",
    file = "data/seqfish.xlsx"
  )
)

# ── Loop over datasets ────────────────────────────────────────────────────────
for (dataset_name in names(datasets)) {
  info <- datasets[[dataset_name]]
  cat(sprintf("\n%s\nProcessing dataset: %s\n%s\n",
              strrep("=", 50), dataset_name, strrep("=", 50)))

  # Check if both scaled and unscaled results already exist
  file_unscaled <- file.path(output_dir, sprintf("SPARK_%s_unscaled_results.csv", dataset_name))
  file_scaled   <- file.path(output_dir, sprintf("SPARK_%s_scaled_results.csv",   dataset_name))

  if (file.exists(file_unscaled) && file.exists(file_scaled)) {
    cat(sprintf("Skipping %s: Both scaled and unscaled results already exist.\n", dataset_name))
    next
  }

  # ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
  if (info$type == "tsv") {
    if (!file.exists(info$expr_file) || !file.exists(info$coord_file)) {
      cat(sprintf("Skipping %s: Files not found.\n", dataset_name))
      next
    }

    # Load expression
    df <- read.delim(info$expr_file, sep = "\t", row.names = "gene", check.names = FALSE)
    if ("ensemblid" %in% colnames(df)) {
      df <- df[, colnames(df) != "ensemblid"]
    }
    counts <- t(df)  # spots x genes

    # Load coordinates
    coord_df    <- read.delim(info$coord_file, sep = "\t", row.names = 1, check.names = FALSE)
    sample_info <- as.data.frame(t(coord_df))

  } else if (info$type == "xlsx") {
    if (!file.exists(info$file)) {
      cat(sprintf("Skipping %s: File not found.\n", dataset_name))
      next
    }

    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    library(readxl)

    # Load expression
    counts_df <- as.data.frame(read_excel(info$file, sheet = "Hippocampus Counts", col_names = FALSE))
    rownames(counts_df) <- counts_df[[1]]
    counts_df <- counts_df[, -1]
    counts <- t(counts_df)  # spots x genes

    # Load coordinates
    coord_df    <- as.data.frame(read_excel(info$file, sheet = "Centroids", col_names = FALSE))
    sample_info <- coord_df[, 1:2]
    colnames(sample_info) <- c("x", "y")
    rownames(sample_info) <- rownames(counts)

  } else {
    cat(sprintf("Unknown file type for %s. Skipping.\n", dataset_name))
    next
  }

  # ── 2. FILTER ────────────────────────────────────────────────────────────────
  # Filter lowly expressed genes (sum >= 3)
  counts <- counts[, colSums(counts) >= 3]

  # Calculate total counts per spot (sequencing depth)
  sample_info$x            <- as.numeric(sample_info$x)
  sample_info$y            <- as.numeric(sample_info$y)
  sample_info$total_counts <- rowSums(counts)

  # Filter out poor quality spots
  valid_spots <- sample_info$total_counts > 10
  sample_info <- sample_info[valid_spots, ]
  counts      <- counts[valid_spots, ]

  # Align coordinates with expression matrix
  sample_info <- sample_info[rownames(counts), ]

  cat(sprintf("Filtered matrix: %d spots, %d genes.\n", nrow(counts), ncol(counts)))

  # ── 3. RUN SPARK (SCALED VS UNSCALED) ────────────────────────────────────────
  for (is_scaled in c(FALSE, TRUE)) {
    scale_label     <- ifelse(is_scaled, "scaled", "unscaled")
    output_filename <- file.path(output_dir,
                                 sprintf("SPARK_%s_%s_results.csv", dataset_name, scale_label))

    if (file.exists(output_filename)) {
      cat(sprintf("\n-> Skipping SPARK (%s): File already exists.\n", scale_label))
      next
    }

    cat(sprintf("\n--- Running SPARK for %s (%s) ---\n", dataset_name, scale_label))

    # Copy coordinates so we don't permanently alter them
    current_sample_info <- sample_info

    if (is_scaled) {
      cat("Scaling coordinates to [0, 1] range...\n")
      current_sample_info$x <- (current_sample_info$x - min(current_sample_info$x)) /
                               (max(current_sample_info$x) - min(current_sample_info$x))
      current_sample_info$y <- (current_sample_info$y - min(current_sample_info$y)) /
                               (max(current_sample_info$y) - min(current_sample_info$y))
    } else {
      cat("Using raw coordinates...\n")
    }

    # SPARK expects genes x spots
    counts_spark <- t(counts)

    # ── Create SPARK object ─────────────────────────────────────────────────
    spark <- CreateSPARKObject(
      counts           = counts_spark,
      location         = current_sample_info[, c("x", "y")],
      percentage       = 0.1,
      min_total_counts = 3
    )
    spark@lib_size <- apply(spark@counts, 2, sum)

    # ── Fit model and test ──────────────────────────────────────────────────
    cat("Fitting variance components...\n")
    spark <- spark.vc(
      spark,
      covariates = NULL,
      lib_size   = spark@lib_size,
      num_core   = 4,
      verbose    = FALSE
    )

    cat("Running hypothesis tests...\n")
    spark <- spark.test(
      spark,
      check_positive = TRUE,
      verbose        = FALSE
    )

    # ── Save results ────────────────────────────────────────────────────────
    results         <- spark@res_mtest
    results$gene    <- rownames(results)
    results_sorted  <- results[order(results$adjusted_pvalue), ]

    cat(sprintf("\nTop 5 Spatially Variable Genes (%s):\n", scale_label))
    print(head(results_sorted[, c("gene", "combined_pvalue", "adjusted_pvalue")], 5))

    write.csv(results_sorted, output_filename, row.names = FALSE)
    cat(sprintf("Saved results to %s\n", output_filename))
  }
}

cat("\nAll datasets processed.\n")