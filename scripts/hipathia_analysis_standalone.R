#!/usr/bin/env Rscript
# HiPathia pathway analysis: all samples vs mock
# Standalone script (no Snakemake). Uses human gene counts and Samplesheet.
# Run from project root: Rscript scripts/hipathia_analysis_standalone.R
# Or: conda activate hipathia && Rscript scripts/hipathia_analysis_standalone.R

# --- Paths ---
# Detect project root: prefer working directory, then look for known marker files
# in parent directories. No hardcoded absolute paths.
find_project_root <- function() {
  if (file.exists("Snakefile") || file.exists("config.yaml") || file.exists("counts/gene_counts_matrix_human.tsv")) {
    return(getwd())
  }
  d <- getwd()
  for (i in 1:5) {
    d <- dirname(d)
    if (file.exists(file.path(d, "Snakefile"))) return(d)
  }
  return(getwd())  # fallback
}
project_dir <- find_project_root()

counts_file  <- file.path(project_dir, "counts", "gene_counts_matrix_human.tsv")
samples_file <- file.path(project_dir, "Samplesheet.tsv")
out_dir      <- file.path(project_dir, "results", "hipathia")
pathway_tsv  <- file.path(out_dir, "hipathia_pathway_activity.tsv")
heatmap_file <- file.path(out_dir, "hipathia_pathway_heatmap.png")
pca_file     <- file.path(out_dir, "hipathia_pca.png")
report_dir   <- file.path(out_dir, "report")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(hipathia)
library(ggplot2)

message("Loading counts and metadata...")
# Counts: gene_id column + sample columns (sample IDs = Name in Samplesheet)
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.table(samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Samplesheet has Sample and Name; count matrix columns are Name (e.g. Unknown_CR065-001T0001)
if (!"Name" %in% colnames(metadata))
  stop("Samplesheet must have a 'Name' column matching count matrix column names.")
rownames(metadata) <- metadata$Name

# Define group: mock vs all others (treated)
metadata$Group <- ifelse(grepl("mock", metadata$Sample, ignore.case = TRUE), "mock", "treated")

# Align to common samples
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0)
  stop("No common samples between counts and Samplesheet (check Name column).")

counts   <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]
message(paste("Using", ncol(counts), "samples and", nrow(counts), "genes"))

# Remove Ensembl version for HiPathia/annotation compatibility
rownames(counts) <- sub("\\.[0-9]+$", "", rownames(counts))
expr_matrix <- as.matrix(counts)

message("Translating gene IDs (Ensembl -> Entrez) for HiPathia...")
trans_data <- translate_data(expr_matrix, species = "hsa", verbose = TRUE)

message("Normalizing expression for HiPathia [0,1]...")
# Use rank-based normalization (by_quantiles=FALSE) to avoid preprocessCore threading issues in some conda envs
exp_data <- normalize_data(trans_data, by_quantiles = FALSE, truncation_percentil = 0.95)

message("Loading pathway data (Homo sapiens)...")
pathways <- load_pathways("hsa")

message("Computing pathway activity...")
pathway_activity <- hipathia(exp_data, pathways, decompose = FALSE, verbose = TRUE)

subpathway_values <- get_paths_data(pathway_activity)
message(paste("Subpathways:", nrow(subpathway_values)))

# Normalize path values
subpathway_norm <- normalize_paths(subpathway_values, pathways)
if (inherits(subpathway_norm, "SummarizedExperiment")) {
  subpathway_values_norm <- assay(subpathway_norm)
} else {
  subpathway_values_norm <- as.matrix(subpathway_norm)
}

# Aggregate subpathways to pathway level
pathway_list <- get_pathways_list(pathways)
subpathway_to_pathway <- sub("^P-(hsa[0-9]+)-.*$", "\\1", rownames(subpathway_values_norm))
if (any(subpathway_to_pathway == rownames(subpathway_values_norm)))
  subpathway_to_pathway <- sub("^.*(hsa[0-9]+).*$", "\\1", rownames(subpathway_values_norm))

pathway_values_list <- list()
for (pid in pathway_list) {
  idx <- which(subpathway_to_pathway == pid)
  if (length(idx) > 0) {
    pathway_values_list[[pid]] <- colMeans(subpathway_values_norm[rownames(subpathway_values_norm)[idx], , drop = FALSE])
  }
}

pathway_values <- do.call(rbind, pathway_values_list)
rownames(pathway_values) <- names(pathway_values_list)
message(paste("Aggregated to", nrow(pathway_values), "pathways"))

# Pathway table with names
pathway_names <- sapply(rownames(pathway_values), function(pid) {
  tryCatch({
    n <- get_path_names(pathways, pid)
    if (length(n) > 0 && !is.na(n) && n != "") n else pid
  }, error = function(e) pid)
})

pathway_df <- as.data.frame(pathway_values)
pathway_df$Pathway <- rownames(pathway_df)
pathway_df$PathwayName <- pathway_names
pathway_df <- pathway_df[, c("Pathway", "PathwayName", setdiff(colnames(pathway_df), c("Pathway", "PathwayName")))]
write.table(pathway_df, pathway_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved:", pathway_tsv))

# Group vector aligned to matrix columns (same order as subpathway_values_norm)
sample_groups <- metadata$Group[match(colnames(subpathway_values_norm), rownames(metadata))]
message(paste("Groups:", paste(unique(sample_groups), collapse = ", ")))

# Per-pathway Wilcoxon test (treated vs mock) with BH correction, using pathway-level activity
message("Running per-pathway Wilcoxon tests (treated vs mock)...")
stats_list <- lapply(seq_len(nrow(pathway_values)), function(i) {
  pid <- rownames(pathway_values)[i]
  pname <- pathway_names[i]
  vals <- pathway_values[i, ]
  g1 <- vals[sample_groups == "treated"]
  g2 <- vals[sample_groups == "mock"]
  if (length(g1) < 2 || length(g2) < 2) {
    return(data.frame(
      Pathway = pid,
      PathwayName = pname,
      n_treated = length(g1),
      n_mock = length(g2),
      median_treated = ifelse(length(g1) > 0, median(g1), NA_real_),
      median_mock = ifelse(length(g2) > 0, median(g2), NA_real_),
      diff_median = NA_real_,
      W = NA_real_,
      p_value = NA_real_
    ))
  }
  wt <- suppressWarnings(wilcox.test(g1, g2, exact = FALSE))
  data.frame(
    Pathway = pid,
    PathwayName = pname,
    n_treated = length(g1),
    n_mock = length(g2),
    median_treated = median(g1),
    median_mock = median(g2),
    diff_median = median(g1) - median(g2),
    W = unname(wt$statistic),
    p_value = wt$p.value
  )
})
pathway_stats <- do.call(rbind, stats_list)
pathway_stats$padj_BH <- p.adjust(pathway_stats$p_value, method = "BH")
stats_tsv <- file.path(out_dir, "hipathia_pathway_stats.tsv")
write.table(pathway_stats, stats_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved pathway stats with BH-adjusted p-values:", stats_tsv))

# Wilcoxon: treated vs mock
comp <- do_wilcoxon(subpathway_values_norm, sample_groups, g1 = "treated", g2 = "mock")

# PCA on pathway activity (prcomp handles more variables than samples; do_pca/princomp does not)
message("Computing PCA on pathway activity...")
subpathway_values_norm[!is.finite(subpathway_values_norm)] <- 0
pathway_mat <- pathway_values  # 146 pathways x 12 samples
pathway_mat[!is.finite(pathway_mat)] <- 0
# Remove constant/zero-variance pathways so scaling works
pathway_sd <- apply(pathway_mat, 1, sd)
pathway_mat <- pathway_mat[pathway_sd > 1e-10, , drop = FALSE]
pca_fit <- prcomp(t(pathway_mat), center = TRUE, scale. = TRUE)  # t() -> samples x variables
pc_scores <- pca_fit$x
var_pct <- 100 * summary(pca_fit)$importance["Proportion of Variance", ]
png(pca_file, width = 800, height = 700, res = 150)
par(mar = c(4, 4, 3, 1))
plot(pc_scores[, 1], pc_scores[, 2], pch = 19, cex = 2, col = as.numeric(factor(sample_groups)),
     xlab = paste0("PC1 (", round(var_pct[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_pct[2], 1), "%)"),
     main = "HiPathia pathway activity: PCA")
legend("bottomleft", legend = unique(sample_groups), col = seq_along(unique(sample_groups)), pch = 19, bty = "n")
text(pc_scores[, 1], pc_scores[, 2], labels = colnames(pathway_mat), pos = 4, cex = 0.7)
dev.off()
message(paste("Saved:", pca_file))

# Heatmap (data already cleaned for PCA)
png(heatmap_file, width = 2000, height = min(3000, 200 + nrow(subpathway_values_norm) * 2), res = 150)
tryCatch({
  heatmap_plot(subpathway_values_norm,
               group = sample_groups,
               colors = "hipathia",
               sample_clust = TRUE,
               variable_clust = FALSE,
               scale = TRUE,
               main = "HiPathia: treated vs mock")
}, error = function(e) {
  heatmap_plot(subpathway_values_norm,
              group = sample_groups,
              colors = "hipathia",
              sample_clust = FALSE,
              variable_clust = FALSE,
              scale = TRUE,
              main = "HiPathia: treated vs mock")
})
dev.off()
message(paste("Saved:", heatmap_file))

# Interactive HTML report (pathway activities visualised on pathway maps in the browser)
# Requires HiPathia's "pathway-viewer" files; conda install may not include them.
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_ok <- FALSE
tryCatch({
  create_report(comp, pathways, output_folder = basename(report_dir), path = dirname(report_dir), verbose = TRUE)
  report_ok <- TRUE
  message(paste("Report saved in:", report_dir))
  message("  Open index.html in that folder in a web browser to view pathway activities.")
  message("  Or in R: hipathia::visualize_report('", report_dir, "')", sep = "")
}, error = function(e) {
  # Try with output_folder = full path (some versions use folder path directly)
  tryCatch({
    create_report(comp, pathways, output_folder = report_dir, verbose = TRUE)
    report_ok <- TRUE
    message(paste("Report saved in:", report_dir))
    message("  Open index.html in a web browser, or in R: hipathia::visualize_report('", report_dir, "')", sep = "")
  }, error = function(e2) {
    message("Interactive report skipped: ", e2$message)
    message("  To get the HTML report, install HiPathia from Bioconductor: BiocManager::install('hipathia')")
  })
})
if (!report_ok)
  message("Results are still in: pathway TSV, PCA plot and heatmap PNG.")

message("HiPathia analysis done.")
