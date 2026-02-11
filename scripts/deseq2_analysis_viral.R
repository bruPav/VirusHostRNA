# DESeq2 Analysis for Viral Genes Only
# This script is similar to deseq2_analysis.R but adapted for viral genes
# which don't have Ensembl IDs, so gene symbol mapping is skipped

# Redirect output to a log file if Snakemake is capturing it
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(vsn)

# 1. Load Data
# Handle empty viral count matrix (0 genes detected)
count_lines <- readLines(snakemake@input$counts)
# Check if file has only a header or is completely empty
if (length(count_lines) <= 1) {
    message("Viral count matrix is empty (0 genes). Creating placeholder outputs.")
    # Write empty results files
    empty_df <- data.frame(GeneID=character(), GeneSymbol=character(), GeneBiotype=character(),
                           baseMean=numeric(), log2FoldChange=numeric(), lfcSE=numeric(),
                           stat=numeric(), pvalue=numeric(), padj=numeric())
    write.table(empty_df, snakemake@output$results, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(empty_df, snakemake@output$results_filtered, sep="\t", quote=FALSE, row.names=FALSE)
    # Write empty normalized counts (just header from count matrix)
    writeLines(count_lines[1], snakemake@output$normalized_counts)
    # Create placeholder plots
    for (plot_out in c(snakemake@output$pca, snakemake@output$distance,
                       snakemake@output$volcano, snakemake@output$heatmap)) {
        png(plot_out, width=800, height=600, res=150)
        plot.new(); title("No viral genes detected in count matrix")
        dev.off()
    }
    message("Placeholder outputs created. Pipeline can continue.")
    sink(type="message"); sink()
    quit(save="no", status=0)
}

counts <- read.table(snakemake@input$counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.table(snakemake@input$samples, header=TRUE, row.names=1, sep="\t")

# Find common samples between counts and metadata
common_samples <- intersect(colnames(counts), rownames(metadata))

if (length(common_samples) == 0) {
    stop("No common samples found between counts matrix and metadata file!")
}

if (length(common_samples) < nrow(metadata)) {
    warning(paste("Only", length(common_samples), "out of", nrow(metadata), "samples from metadata are present in counts matrix."))
    warning(paste("Missing samples:", paste(setdiff(rownames(metadata), common_samples), collapse=", ")))
}

# Filter both to only include common samples and ensure same order
counts <- counts[, common_samples, drop=FALSE]
metadata <- metadata[common_samples, , drop=FALSE]

# Check if we have enough genes (viral genes are few, so we need at least 2)
if (nrow(counts) < 2) {
    warning("Very few viral genes detected. DESeq2 analysis may be limited.")
}

# 2. Create DESeq2 Dataset with simple design
# Use Control column if available, otherwise use Condition
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
    stop("Metadata must contain either 'Control' or 'Condition' column.")
}

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = design_formula)

# Viral counts are very sparse (many zeros); default size factor estimation fails.
# Compute size factors with poscounts, then replace any NA (e.g. zero viral samples) before assigning.
cnts <- counts(dds)
sf <- DESeq2::estimateSizeFactorsForMatrix(cnts, type = "poscounts")
if (any(is.na(sf))) {
  message("Some size factors were NA (e.g. zero viral counts); setting them to 1.")
  sf[is.na(sf)] <- 1
}
sizeFactors(dds) <- sf
# For small gene sets, use betaPrior=FALSE to avoid convergence issues
dds <- DESeq(dds, betaPrior=FALSE)

# 3. Extract Results
if ("Control" %in% colnames(metadata) && "Experiment" %in% levels(dds$Control)) {
    res <- results(dds, contrast=c("Control", "Experiment", "Negative"))
} else {
    res <- results(dds)
}

# Apply Benjamini-Hochberg is done by default in results()
res_df <- as.data.frame(res)
# Convert row names to a column with proper header
res_df$GeneID <- rownames(res_df)

# For viral genes, we use the gene ID as the symbol (no Ensembl mapping needed)
res_df$GeneSymbol <- res_df$GeneID
res_df$GeneBiotype <- "viral_gene"

# Reorder columns: GeneID, GeneSymbol, GeneBiotype, then statistics
res_df <- res_df[, c("GeneID", "GeneSymbol", "GeneBiotype", setdiff(colnames(res_df), c("GeneID", "GeneSymbol", "GeneBiotype")))]
write.table(res_df, snakemake@output$results, sep="\t", quote=FALSE, row.names=FALSE)

# Filter and write significant genes (padj < 0.05)
res_df_significant <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
message(paste("Found", nrow(res_df_significant), "significantly differentially expressed viral genes (padj < 0.05)"))
write.table(res_df_significant, snakemake@output$results_filtered, sep="\t", quote=FALSE, row.names=FALSE)

# Save normalized counts for downstream analysis
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$GeneID <- rownames(normalized_counts_df)
normalized_counts_df <- normalized_counts_df[, c("GeneID", setdiff(colnames(normalized_counts_df), "GeneID"))]
write.table(normalized_counts_df, snakemake@output$normalized_counts, sep="\t", quote=FALSE, row.names=FALSE)
message("Saved normalized counts for viral genes")

# 4. Transformations for Visualization
# With few sparse viral genes, vst() often fails (needs enough genes with mean count > 5).
# Use vst if possible, else fall back to log2(normalized + 1) for plots.
n_genes <- nrow(dds)
vsd <- tryCatch(vst(dds, blind = FALSE, nsub = min(1000L, n_genes)), error = function(e) NULL)
use_vst <- !is.null(vsd)
if (!use_vst) {
  message("VST failed (typical for few viral genes); using log2(normalized counts + 1) for visualizations.")
  norm_counts <- counts(dds, normalized = TRUE)
  transform_mat <- log2(norm_counts + 1)
}

# Build a vector of condition labels aligned to sample columns
sample_names <- if (use_vst) colnames(vsd) else colnames(transform_mat)
condition_labels <- as.character(metadata$Condition[match(sample_names, rownames(metadata))])

# --- VISUALIZATION 1: PCA Plot ---
png(snakemake@output$pca, width=1200, height=1000, res=150)
if (use_vst) {
  pcaData <- plotPCA(vsd, intgroup = c("Condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
} else {
  pca <- prcomp(t(transform_mat), center = TRUE, scale. = FALSE)
  percentVar <- round(100 * summary(pca)$importance[2, 1:2])
  pcaData <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition_labels)
}
ggplot(pcaData, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA: Sample Clustering (Viral Genes)")
dev.off()

# --- VISUALIZATION 2: Sample Distance Matrix ---
if (use_vst) {
  sampleDists <- dist(t(assay(vsd)))
} else {
  sampleDists <- dist(t(transform_mat))
}
sampleDistMatrix <- as.matrix(sampleDists)
display_labels <- paste0(condition_labels, "_", seq_along(condition_labels))
rownames(sampleDistMatrix) <- display_labels
colnames(sampleDistMatrix) <- display_labels
png(snakemake@output$distance, width = 1200, height = 1000, res = 150)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distances (Viral Genes)")
dev.off()

# --- VISUALIZATION 3: Volcano Plot ---
png(snakemake@output$volcano, width=1200, height=1200, res=150)
# Use gene IDs for labels (viral genes don't have symbols)
res_df$gene_label <- res_df$GeneID
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color = padj < 0.05), alpha=0.6, size=3) +
    theme_minimal() +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(data=subset(res_df, !is.na(padj) & padj < 0.05), aes(label=gene_label), size=4) +
    labs(title="Differential Expression Volcano Plot (Viral Genes)", 
         subtitle="Red: Adjusted P-Value < 0.05",
         x="log2 Fold Change",
         y="-log10 Adjusted P-value")
dev.off()

# --- VISUALIZATION 4: Significant DE Genes Heatmap ---
# Use all genes with padj < 0.05
significant_genes_idx <- which(!is.na(res$padj) & res$padj < 0.05)

# Adjust height based on number of genes (minimum 800, add 30 pixels per gene)
n_sig_genes <- length(significant_genes_idx)
heatmap_height <- max(800, 200 + n_sig_genes * 30)

# Use same transform as for PCA/distance (vst or log2)
heatmap_transform <- if (use_vst) assay(vsd) else transform_mat

png(snakemake@output$heatmap, width = 1200, height = heatmap_height, res = 150)
if (n_sig_genes == 0) {
    message("No viral genes with padj < 0.05; skipping heatmap.")
    plot.new()
    title("No significantly differentially expressed viral genes (padj < 0.05)")
} else {
    heatmap_mat <- heatmap_transform[significant_genes_idx, , drop = FALSE]
    # Use gene IDs for row names
    rownames(heatmap_mat) <- rownames(heatmap_mat)
    # Build unique display labels that still show Condition
    display_labels <- paste0(condition_labels, "_", seq_along(condition_labels))
    colnames(heatmap_mat) <- display_labels
    # Build annotation data frame whose rownames match the heatmap columns
    ann_col <- data.frame(Condition = condition_labels)
    rownames(ann_col) <- display_labels
    # hclust needs at least 2 objects; disable clustering when only 1 row or 1 column
    cluster_rows <- n_sig_genes >= 2
    cluster_cols <- ncol(heatmap_mat) >= 2
    pheatmap(heatmap_mat,
             cluster_rows = cluster_rows,
             show_rownames = TRUE,
             cluster_cols = cluster_cols,
             annotation_col = ann_col,
             fontsize_row = 10,
             main = paste("Significantly Differentially Expressed Viral Genes (padj < 0.05, n =", n_sig_genes, ")"))
}
dev.off()

message("Viral gene DESeq2 analysis completed successfully")
