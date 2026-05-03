library(SPARK)
library(Matrix)

output_dir <- "SPARK_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- list(
  Rep11_MOB = list(type = "tsv_counts", file = "data/Rep11_MOB_count_matrix-1.tsv"),
  Layer_2_BC = list(type = "rds",        file = "data/Layer2_BC_Count.rds"),
  seqfish    = list(type = "xlsx",       file = "data/seqfish.xlsx")
)

# Helper: coerce any matrix to a raw integer count matrix.
# Warns if values look non-integer (e.g. log-normalised).
as_raw_counts <- function(mat, name) {
  if (is.character(mat)) {
    mat <- apply(mat, 2, as.numeric)
  }
  storage.mode(mat) <- "double"

  if (any(mat < 0, na.rm = TRUE))
    stop(sprintf("[%s] Negative values detected – matrix does not look like raw counts.", name))

  frac <- mean(abs(mat - round(mat)) > 1e-6, na.rm = TRUE)
  if (frac > 0.01)
    warning(sprintf(
      "[%s] %.1f%% of expression values are non-integer. SPARK requires raw counts — log-normalised or CPM input will produce invalid results.",
      name, frac * 100))

  storage.mode(mat) <- "integer"
  mat
}

for (dataset_name in names(datasets)) {
  info <- datasets[[dataset_name]]
  cat(sprintf("\n%s\nProcessing dataset: %s\n%s\n",
              strrep("=", 50), dataset_name, strrep("=", 50)))

  file_unscaled <- file.path(output_dir, sprintf("SPARK_%s_unscaled_results.csv", dataset_name))
  file_scaled   <- file.path(output_dir, sprintf("SPARK_%s_scaled_results.csv",   dataset_name))

  if (file.exists(file_unscaled) && file.exists(file_scaled)) {
    cat(sprintf("Skipping %s: both result files already exist.\n", dataset_name))
    next
  }

  # ── 1. LOAD ───────────────────────────────────────────────────────────────
  if (info$type == "tsv_counts") {
    if (!file.exists(info$file)) { cat("Skipping: file not found.\n"); next }

    # Read without row.names=1 — first column header is empty in this file,
    # which causes R to silently replace spot IDs with sequential integers.
    df       <- read.delim(info$file, sep = "\t", check.names = FALSE)
    spot_ids <- as.character(df[[1]])
    df       <- df[, -1, drop = FALSE]
    rownames(df) <- spot_ids

    cat(sprintf("  Raw TSV: %d spots × %d genes.\n", nrow(df), ncol(df)))
    cat(sprintf("  Example spot IDs: %s\n", paste(head(spot_ids, 4), collapse = ", ")))

    # Guard: every column must be numeric (columns = genes after the ID column is removed)
    non_num <- names(df)[!sapply(df, function(x) is.numeric(x) || is.integer(x))]
    if (length(non_num))
      stop(sprintf("[%s] Non-numeric gene columns: %s",
                   dataset_name, paste(head(non_num, 10), collapse = ", ")))

    # Transpose to spots × genes, coerce to integer
    counts <- as_raw_counts(as.matrix(df), dataset_name)

    # Parse coordinates from spot IDs: "16.92x9.015" → x=16.92, y=9.015
    spot_names <- rownames(counts)

    bad_splits <- which(sapply(strsplit(spot_names, "x"), length) != 2)
    if (length(bad_splits))
      stop(sprintf("[%s] %d spot names did not split into 2 coordinate parts: %s",
                   dataset_name, length(bad_splits),
                   paste(head(spot_names[bad_splits], 5), collapse = ", ")))

    coords   <- do.call(rbind, strsplit(spot_names, "x"))
    x_coords <- suppressWarnings(as.numeric(coords[, 1]))
    y_coords <- suppressWarnings(as.numeric(coords[, 2]))

    if (any(is.na(x_coords)) || any(is.na(y_coords)))
      stop(sprintf("[%s] NA coordinates after parsing spot names. Example bad names: %s",
                   dataset_name,
                   paste(head(spot_names[is.na(x_coords) | is.na(y_coords)], 5), collapse = ", ")))

    sample_info <- data.frame(x = x_coords, y = y_coords, row.names = spot_names)
    cat(sprintf("  Coordinate ranges: x=[%.3f, %.3f], y=[%.3f, %.3f]\n",
                min(x_coords), max(x_coords), min(y_coords), max(y_coords)))

  } else if (info$type == "rds") {
    if (!file.exists(info$file)) { cat("Skipping: file not found.\n"); next }

    env        <- new.env()
    loaded_obj <- tryCatch(readRDS(info$file), error = function(e) {
      load(info$file, envir = env)
      obj_names <- ls(env)
      if (!length(obj_names)) stop("No objects found in loaded file.")
      env[[obj_names[1]]]
    })

    counts <- as_raw_counts(t(as.matrix(loaded_obj)), dataset_name)

    spot_names <- rownames(counts)

    bad_splits <- which(sapply(strsplit(spot_names, "x"), length) != 2)
    if (length(bad_splits))
      stop(sprintf("[%s] %d spot names did not split into 2 coordinate parts: %s",
                   dataset_name, length(bad_splits),
                   paste(head(spot_names[bad_splits], 5), collapse = ", ")))

    coords   <- do.call(rbind, strsplit(spot_names, "x"))
    x_coords <- suppressWarnings(as.numeric(coords[, 1]))
    y_coords <- suppressWarnings(as.numeric(coords[, 2]))

    if (any(is.na(x_coords)) || any(is.na(y_coords)))
      stop(sprintf("[%s] NA coordinates after parsing. Example bad names: %s",
                   dataset_name,
                   paste(head(spot_names[is.na(x_coords) | is.na(y_coords)], 5), collapse = ", ")))

    sample_info <- data.frame(x = x_coords, y = y_coords, row.names = spot_names)

  } else if (info$type == "xlsx") {
    if (!file.exists(info$file)) { cat("Skipping: file not found.\n"); next }
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    library(readxl)

    # col_names = TRUE: first row is gene names, first column is cell/spot IDs
    counts_raw           <- as.data.frame(readxl::read_excel(info$file, sheet = "Hippocampus Counts", col_names = TRUE))
    rownames(counts_raw) <- counts_raw[[1]]
    counts_raw           <- counts_raw[, -1, drop = FALSE]

    counts_raw[] <- lapply(counts_raw, function(x) suppressWarnings(as.numeric(as.character(x))))
    if (anyNA(counts_raw))
      warning(sprintf("[%s] NAs introduced when coercing count sheet to numeric. Check for non-numeric cells.", dataset_name))

    counts <- as_raw_counts(t(as.matrix(counts_raw)), dataset_name)

    # Coordinates: pick the first two numeric columns
    coord_df     <- as.data.frame(readxl::read_excel(info$file, sheet = "Centroids", col_names = TRUE))
    numeric_cols <- names(coord_df)[sapply(coord_df, is.numeric)]
    if (length(numeric_cols) < 2)
      stop(sprintf("[%s] Centroids sheet has fewer than 2 numeric columns.", dataset_name))

    sample_info <- data.frame(
      x = as.numeric(coord_df[[numeric_cols[1]]]),
      y = as.numeric(coord_df[[numeric_cols[2]]]),
      row.names = rownames(counts)
    )

    if (nrow(sample_info) != nrow(counts))
      stop(sprintf("[%s] Mismatch: %d spots in count matrix but %d rows in Centroids sheet.",
                   dataset_name, nrow(counts), nrow(sample_info)))

  } else {
    cat(sprintf("Unknown file type for %s. Skipping.\n", dataset_name)); next
  }

  # ── 2. FILTER ─────────────────────────────────────────────────────────────
  counts <- counts[, colSums(counts) >= 3L, drop = FALSE]

  sample_info$x            <- as.numeric(sample_info$x)
  sample_info$y            <- as.numeric(sample_info$y)
  sample_info$total_counts <- rowSums(counts)

  valid_spots <- sample_info$total_counts > 10L
  sample_info <- sample_info[valid_spots, , drop = FALSE]
  counts      <- counts[valid_spots, , drop = FALSE]
  sample_info <- sample_info[rownames(counts), , drop = FALSE]

  cat(sprintf("Filtered matrix: %d spots × %d genes.\n", nrow(counts), ncol(counts)))

  # ── 3. SPARK (coordinate-scaled vs unscaled) ───────────────────────────────
  for (is_scaled in c(FALSE, TRUE)) {
    scale_label     <- ifelse(is_scaled, "scaled", "unscaled")
    output_filename <- file.path(output_dir,
                                 sprintf("SPARK_%s_%s_results.csv", dataset_name, scale_label))

    if (file.exists(output_filename)) {
      cat(sprintf("Skipping SPARK (%s): file already exists.\n", scale_label)); next
    }

    cat(sprintf("\n--- Running SPARK for %s (%s coordinates) ---\n", dataset_name, scale_label))
    current_sample_info <- sample_info

    if (is_scaled) {
      cat("Scaling spatial coordinates to [0, 1]...\n")
      rng_x <- range(current_sample_info$x)
      rng_y <- range(current_sample_info$y)
      if (diff(rng_x) == 0 || diff(rng_y) == 0)
        stop(sprintf("[%s] All coordinates identical on one axis – cannot scale.", dataset_name))
      current_sample_info$x <- (current_sample_info$x - rng_x[1]) / diff(rng_x)
      current_sample_info$y <- (current_sample_info$y - rng_y[1]) / diff(rng_y)
    } else {
      cat("Using raw spatial coordinates...\n")
    }

    # SPARK expects genes × spots
    counts_spark <- t(counts)

    spark <- CreateSPARKObject(
      counts           = counts_spark,
      location         = current_sample_info[, c("x", "y")],
      percentage       = 0.1,
      min_total_counts = 3
    )
    spark@lib_size <- colSums(spark@counts)

    cat("Fitting variance components...\n")
    spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size,
                      num_core = 4, verbose = FALSE)

    cat("Running hypothesis tests...\n")
    spark <- spark.test(spark, check_positive = TRUE, verbose = FALSE)

    results        <- spark@res_mtest
    results$gene   <- rownames(results)
    results_sorted <- results[order(results$adjusted_pvalue), ]

    cat(sprintf("\nTop 5 SVGs (%s coordinates):\n", scale_label))
    print(head(results_sorted[, c("gene", "combined_pvalue", "adjusted_pvalue")], 5))

    write.csv(results_sorted, output_filename, row.names = FALSE)
    cat(sprintf("Saved: %s\n", output_filename))
  }
}

cat("\nAll datasets processed.\n")