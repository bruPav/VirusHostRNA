# Redirect output to a log file if Snakemake is capturing it
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(vsn)
library(biomaRt)

# 1. Load Data
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

dds <- DESeq(dds)

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

# Map Ensembl IDs to gene symbols using biomaRt
message("Mapping Ensembl IDs to gene symbols...")
# Remove version numbers from Ensembl IDs (e.g., ENSG00000123456.1 -> ENSG00000123456)
ensembl_ids <- sub("\\.[0-9]+$", "", res_df$GeneID)

# Initialize gene symbol columns with NA
res_df$GeneSymbol <- NA
res_df$GeneBiotype <- NA

# Try to connect to Ensembl database and get gene symbols
tryCatch({
    # Connect to Ensembl database
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Get gene symbols (process in batches to avoid timeout)
    batch_size <- 10000
    gene_map_list <- list()
    for (i in seq(1, length(ensembl_ids), by = batch_size)) {
        end_idx <- min(i + batch_size - 1, length(ensembl_ids))
        batch_ids <- ensembl_ids[i:end_idx]
        
        gene_map_batch <- getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
            filters = "ensembl_gene_id",
            values = batch_ids,
            mart = ensembl
        )
        gene_map_list[[length(gene_map_list) + 1]] <- gene_map_batch
    }
    gene_map <- do.call(rbind, gene_map_list)
    
    # Merge with results (handle cases where multiple symbols map to same ID or no mapping found)
    # Take first symbol if multiple exist
    gene_map_unique <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
    res_df$GeneSymbol <- gene_map_unique$hgnc_symbol[match(ensembl_ids, gene_map_unique$ensembl_gene_id)]
    res_df$GeneBiotype <- gene_map_unique$gene_biotype[match(ensembl_ids, gene_map_unique$ensembl_gene_id)]
    
    # Replace empty strings with NA
    res_df$GeneSymbol[res_df$GeneSymbol == ""] <- NA
    
    message(paste("Successfully mapped", sum(!is.na(res_df$GeneSymbol)), "genes to symbols"))
}, error = function(e) {
    warning(paste("Failed to map gene symbols using biomaRt:", e$message))
    warning("Results will contain only Ensembl IDs")
})

# Reorder columns: GeneID, GeneSymbol, GeneBiotype, then statistics
res_df <- res_df[, c("GeneID", "GeneSymbol", "GeneBiotype", setdiff(colnames(res_df), c("GeneID", "GeneSymbol", "GeneBiotype")))]
write.table(res_df, snakemake@output$results, sep="\t", quote=FALSE, row.names=FALSE)

# Filter and write significant genes (padj < 0.05)
res_df_significant <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
message(paste("Found", nrow(res_df_significant), "significantly differentially expressed genes (padj < 0.05)"))
write.table(res_df_significant, snakemake@output$results_filtered, sep="\t", quote=FALSE, row.names=FALSE)

# Save normalized counts for downstream analysis (e.g., HiPathia)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$GeneID <- rownames(normalized_counts_df)
normalized_counts_df <- normalized_counts_df[, c("GeneID", setdiff(colnames(normalized_counts_df), "GeneID"))]
write.table(normalized_counts_df, snakemake@output$normalized_counts, sep="\t", quote=FALSE, row.names=FALSE)
message("Saved normalized counts for downstream pathway analysis")

# 4. Transformations for Visualization
vsd <- vst(dds, blind=FALSE)

# Build a vector of condition labels aligned to the columns/samples in 'vsd'
condition_labels <- as.character(metadata$Condition[match(colnames(vsd), rownames(metadata))])

# --- VISUALIZATION 1: PCA Plot ---
png(snakemake@output$pca, width=1200, height=1000, res=150)
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=4, alpha=0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA: Sample Clustering")
dev.off()

# --- VISUALIZATION 2: Sample Distance Matrix ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# Relabel rows and columns by Condition (with replicate index) instead of sample IDs
display_labels <- paste0(condition_labels, "_", seq_along(condition_labels))
rownames(sampleDistMatrix) <- display_labels
colnames(sampleDistMatrix) <- display_labels
png(snakemake@output$distance, width=1200, height=1000, res=150)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample-to-Sample Distances")
dev.off()

# --- VISUALIZATION 3: Volcano Plot ---
png(snakemake@output$volcano, width=1200, height=1200, res=150)
# Use gene symbols for labels, fallback to GeneID if symbol is missing
res_df$gene_label <- ifelse(!is.na(res_df$GeneSymbol), res_df$GeneSymbol, res_df$GeneID)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color = padj < 0.05), alpha=0.4) +
    theme_minimal() +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(data=subset(res_df, padj < 0.000001), aes(label=gene_label), size=3) +
    labs(title="Differential Expression Volcano Plot", 
         subtitle="Red: Adjusted P-Value < 0.05",
         x="log2 Fold Change",
         y="-log10 Adjusted P-value")
dev.off()

# --- VISUALIZATION 4: Significant DE Genes Heatmap ---
# Use all genes with padj < 0.05
significant_genes_idx <- which(!is.na(res$padj) & res$padj < 0.05)

# Adjust height based on number of genes (minimum 800, add 20 pixels per gene)
n_sig_genes <- length(significant_genes_idx)
heatmap_height <- max(800, 200 + n_sig_genes * 20)

png(snakemake@output$heatmap, width=1200, height=heatmap_height, res=150)
if (n_sig_genes == 0) {
    message("No genes with padj < 0.05; skipping heatmap.")
    plot.new()
    title("No significantly differentially expressed genes (padj < 0.05)")
} else {
    heatmap_mat <- assay(vsd)[significant_genes_idx, , drop = FALSE]
    # Use gene symbols for row names, fallback to Ensembl ID if symbol is missing
    gene_ids_sig <- rownames(heatmap_mat)
    gene_symbols_sig <- res_df$GeneSymbol[match(gene_ids_sig, res_df$GeneID)]
    gene_labels_sig <- ifelse(!is.na(gene_symbols_sig), gene_symbols_sig, gene_ids_sig)
    rownames(heatmap_mat) <- gene_labels_sig
    # Build unique display labels that still show Condition
    display_labels <- paste0(condition_labels, "_", seq_along(condition_labels))
    colnames(heatmap_mat) <- display_labels
    # Build annotation data frame whose rownames match the heatmap columns
    ann_col <- data.frame(Condition = condition_labels)
    rownames(ann_col) <- display_labels
    pheatmap(heatmap_mat, 
             cluster_rows=TRUE, 
             show_rownames=TRUE,
             cluster_cols=TRUE, 
             annotation_col=ann_col,
             fontsize_row = ifelse(n_sig_genes > 100, 6, 8),
             main=paste("Significantly Differentially Expressed Genes (padj < 0.05, n =", n_sig_genes, ")"))
}
dev.off()
