# =============================================================================
# Shared DESeq2 pipeline core — sourced by deseq2_analysis.R and
# deseq2_analysis_viral.R. Defines reusable functions for log setup, data
# loading, design, pre-filtering, DESeq2 execution, and visualization.
# =============================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# ---------------------------------------------------------------------------
# 1. Log setup
# ---------------------------------------------------------------------------
setup_log <- function(log_path) {
  log <- file(log_path, open = "wt")
  sink(log)
  sink(log, type = "message")
  log
}

# ---------------------------------------------------------------------------
# 2. Load count matrix and metadata, align on common samples
# ---------------------------------------------------------------------------
load_common_samples <- function(counts_file, metadata_file) {
  counts <- read.table(counts_file, header = TRUE, row.names = 1,
                        sep = "\t", check.names = FALSE)
  metadata <- read.table(metadata_file, header = TRUE, row.names = 1,
                          sep = "\t")
  common <- intersect(colnames(counts), rownames(metadata))
  if (length(common) == 0) {
    stop("No common samples found between counts matrix and metadata file!")
  }
  if (length(common) < nrow(metadata)) {
    msg <- paste("Only", length(common), "of", nrow(metadata),
                 "metadata samples present in counts matrix.")
    warning(msg)
    warning(paste("Missing:", paste(setdiff(rownames(metadata), common),
                                    collapse = ", ")))
  }
  list(counts   = counts[, common, drop = FALSE],
       metadata = metadata[common, , drop = FALSE])
}

# ---------------------------------------------------------------------------
# 3. Determine design formula, condition column, and thresholds from config
# ---------------------------------------------------------------------------
setup_design <- function(metadata, config) {
  design_str <- config[["design_formula"]]
  if (!is.null(design_str) && nchar(design_str) > 0) {
    design_formula <- as.formula(design_str)
    for (col in all.vars(design_formula)) {
      if (col %in% colnames(metadata)) {
        metadata[[col]] <- factor(metadata[[col]])
      } else {
        stop(paste("Design formula column", col, "not found in metadata"))
      }
    }
  } else {
    if ("Control" %in% colnames(metadata)) {
      metadata$Control <- factor(metadata$Control)
      if ("Negative" %in% levels(metadata$Control)) {
        metadata$Control <- relevel(metadata$Control, ref = "Negative")
      }
      design_formula <- ~ Control
    } else if ("Condition" %in% colnames(metadata)) {
      metadata$Condition <- factor(metadata$Condition)
      design_formula <- ~ Condition
    } else {
      stop("Metadata needs Control/Condition column, or set design_formula in config")
    }
  }
  condition_col <- if ("Condition" %in% colnames(metadata)) "Condition" else "Control"
  padj_threshold <- if (is.null(config[["deseq2_padj"]])) 0.05 else as.numeric(config[["deseq2_padj"]])
  list(metadata        = metadata,
       design_formula  = design_formula,
       condition_col   = condition_col,
       padj_threshold  = padj_threshold)
}

# ---------------------------------------------------------------------------
# 4. Pre-filter low-count genes
# ---------------------------------------------------------------------------
prefilter_counts <- function(counts, config, default_min_count = 10,
                              default_min_samples = 3) {
  mc <- if (is.null(config[["min_count"]])) default_min_count else as.numeric(config[["min_count"]])
  ms <- if (is.null(config[["min_samples"]])) min(default_min_samples, ncol(counts)) else as.numeric(config[["min_samples"]])
  keep <- rowSums(counts >= mc) >= ms
  message(paste("Pre-filtering: kept", sum(keep), "of", nrow(counts), "genes",
                "(min", mc, "counts in", ms, "samples)"))
  counts[keep, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 5. Build DESeq2 dataset and run DESeq
# ---------------------------------------------------------------------------
run_deseq2 <- function(counts, metadata, design_formula) {
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData   = metadata,
                                design    = design_formula)
  dds <- DESeq(dds)
  dds
}

# ---------------------------------------------------------------------------
# 6. Extract results with optional contrast from config
# ---------------------------------------------------------------------------
extract_results <- function(dds, metadata, config) {
  contrast_cfg <- config[["deseq2_contrast"]]
  if (!is.null(contrast_cfg) && length(contrast_cfg) == 3) {
    results(dds, contrast = contrast_cfg)
  } else if ("Control" %in% colnames(metadata) && "Experiment" %in% levels(dds$Control)) {
    results(dds, contrast = c("Control", "Experiment", "Negative"))
  } else {
    results(dds)
  }
}

# ---------------------------------------------------------------------------
# 7. Write DE results, significant table, and normalized counts
# ---------------------------------------------------------------------------
write_outputs <- function(res_df, dds, padj_threshold, output) {
  res_df$GeneID <- rownames(res_df)
  cols <- c("GeneID", "GeneSymbol", "GeneBiotype",
            setdiff(colnames(res_df), c("GeneID", "GeneSymbol", "GeneBiotype")))
  res_df <- res_df[, cols]
  write.table(res_df, output$results, sep = "\t", quote = FALSE, row.names = FALSE)

  res_sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold, ]
  message(paste("Found", nrow(res_sig),
                "significantly DE genes (padj <", padj_threshold, ")"))
  write.table(res_sig, output$results_filtered, sep = "\t", quote = FALSE, row.names = FALSE)

  nrm <- as.data.frame(counts(dds, normalized = TRUE))
  nrm$GeneID <- rownames(nrm)
  nrm <- nrm[, c("GeneID", setdiff(colnames(nrm), "GeneID"))]
  write.table(nrm, output$normalized_counts, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved normalized counts for downstream analysis")
}

# ---------------------------------------------------------------------------
# 8-11. Visualization helpers
# ---------------------------------------------------------------------------

# Build a vector of condition labels aligned to sample columns
get_condition_labels <- function(transform_data, metadata, condition_col,
                                  use_vst = TRUE) {
  sample_names <- colnames(if (use_vst) assay(transform_data) else transform_data)
  as.character(metadata[[condition_col]][match(sample_names, rownames(metadata))])
}

# PCA plot
plot_pca <- function(transform_data, metadata, condition_col, outfile,
                      title_text = "", use_vst = TRUE) {
  png(outfile, width = 1200, height = 1000, res = 150)
  if (use_vst) {
    pcaData <- plotPCA(transform_data, intgroup = condition_col,
                        returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
  } else {
    mat <- transform_data
    pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
    percentVar <- round(100 * summary(pca)$importance[2, 1:2])
    cl <- as.character(metadata[[condition_col]][match(colnames(mat), rownames(metadata))])
    pcaData <- setNames(
      data.frame(pca$x[, 1], pca$x[, 2], cl, stringsAsFactors = FALSE),
      c("PC1", "PC2", condition_col))
  }
  p <- ggplot(pcaData, aes(PC1, PC2, color = .data[[condition_col]])) +
    geom_point(size = 4, alpha = 0.8) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    ggtitle(paste0("PCA: Sample Clustering", title_text))
  print(p)
  dev.off()
}

# Sample distance matrix heatmap
plot_distance_matrix <- function(transform_data, metadata, condition_col,
                                  outfile, title_text = "", use_vst = TRUE) {
  mat <- if (use_vst) assay(transform_data) else transform_data
  cl  <- as.character(metadata[[condition_col]][match(colnames(mat), rownames(metadata))])
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  display_labels <- paste0(cl, "_", seq_along(cl))
  rownames(sampleDistMatrix) <- display_labels
  colnames(sampleDistMatrix) <- display_labels
  png(outfile, width = 1200, height = 1000, res = 150)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           main = paste0("Sample-to-Sample Distances", title_text))
  dev.off()
}

# Volcano plot
plot_volcano <- function(res_df, padj_threshold, gene_label_col, outfile,
                          title_text = "", label_threshold = 0.000001,
                          point_alpha = 0.4) {
  png(outfile, width = 1200, height = 1200, res = 150)
  res_df$gene_label <- ifelse(!is.na(res_df[[gene_label_col]]),
                                res_df[[gene_label_col]], res_df$GeneID)
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < padj_threshold), alpha = point_alpha) +
    theme_minimal() +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(
      data = subset(res_df, !is.na(padj) & padj < label_threshold),
      aes(label = gene_label), size = 3) +
    labs(title    = paste0("Differential Expression Volcano Plot", title_text),
         subtitle = paste0("Red: Adjusted P-Value < ", padj_threshold),
         x = "log2 Fold Change", y = "-log10 Adjusted P-value")
  print(p)
  dev.off()
}

# Heatmap of significant DE genes
plot_heatmap <- function(transform_data, res, res_df, padj_threshold,
                          metadata, condition_col, gene_label_col, outfile,
                          title_text = "", use_vst = TRUE) {
  sig_idx <- which(!is.na(res$padj) & res$padj < padj_threshold)
  n_sig   <- length(sig_idx)
  hm_height <- max(800, 200 + n_sig * 20)
  hm_mat <- if (use_vst) assay(transform_data) else transform_data
  cl <- as.character(metadata[[condition_col]][match(colnames(hm_mat), rownames(metadata))])

  png(outfile, width = 1200, height = hm_height, res = 150)
  if (n_sig == 0) {
    message(paste0("No genes with padj < ", padj_threshold, "; skipping heatmap."))
    plot.new()
    title(paste0("No significantly DE genes (padj < ", padj_threshold, ")", title_text))
  } else {
    hm <- hm_mat[sig_idx, , drop = FALSE]
    gene_ids <- rownames(hm)
    gene_sym <- res_df[[gene_label_col]][match(gene_ids, res_df$GeneID)]
    gene_lbl <- ifelse(!is.na(gene_sym), gene_sym, gene_ids)
    rownames(hm) <- gene_lbl
    disp_labels <- paste0(cl, "_", seq_along(cl))
    colnames(hm) <- disp_labels
    ann_col <- setNames(data.frame(cl, stringsAsFactors = FALSE), condition_col)
    rownames(ann_col) <- disp_labels
    cr <- n_sig >= 2
    cc <- ncol(hm) >= 2
    pheatmap(hm,
             cluster_rows    = cr,
             show_rownames   = TRUE,
             cluster_cols    = cc,
             annotation_col  = ann_col,
             fontsize_row    = ifelse(n_sig > 100, 6, 8),
             main = paste0("Significantly DE Genes (padj < ", padj_threshold,
                           ", n = ", n_sig, ")", title_text))
  }
  dev.off()
}

# ---------------------------------------------------------------------------
# 12. Run all four visualizations at once
# ---------------------------------------------------------------------------
make_visualizations <- function(transform_data, res_df, res, metadata,
                                 condition_col, padj_threshold, gene_label_col,
                                 output, title_text = "", use_vst = TRUE,
                                 label_threshold = 0.000001,
                                 volcano_alpha = 0.4) {
  cl <- get_condition_labels(transform_data, metadata, condition_col, use_vst)
  plot_pca(transform_data, metadata, condition_col,
            output$pca, title_text, use_vst)
  plot_distance_matrix(transform_data, metadata, condition_col,
                        output$distance, title_text, use_vst)
  plot_volcano(res_df, padj_threshold, gene_label_col,
                output$volcano, title_text, label_threshold, volcano_alpha)
  plot_heatmap(transform_data, res, res_df, padj_threshold,
                metadata, condition_col, gene_label_col,
                output$heatmap, title_text, use_vst)
}

# ---------------------------------------------------------------------------
# 13. Teardown
# ---------------------------------------------------------------------------
teardown_log <- function(log) {
  sink(type = "message")
  sink()
  close(log)
}
