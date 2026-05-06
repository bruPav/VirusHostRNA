# DESeq2 Analysis for Viral Genes Only
# Uses shared core functions from deseq2_core.R with viral-specific
# adaptations: sparse count handling, VST fallback, no Ensembl mapping.
source("scripts/deseq2_core.R")
library(vsn)

# 1. Setup
log_con <- setup_log(snakemake@log[[1]])

# 2. Handle empty viral count matrix (0 genes) — early exit
count_lines <- readLines(snakemake@input$counts)
if (length(count_lines) <= 1) {
  message("Viral count matrix is empty (0 genes). Creating placeholder outputs.")
  empty_df <- data.frame(
    GeneID = character(), GeneSymbol = character(), GeneBiotype = character(),
    baseMean = numeric(), log2FoldChange = numeric(), lfcSE = numeric(),
    stat = numeric(), pvalue = numeric(), padj = numeric())
  write.table(empty_df, snakemake@output$results, sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(empty_df, snakemake@output$results_filtered, sep = "\t",
              quote = FALSE, row.names = FALSE)
  writeLines(count_lines[1], snakemake@output$normalized_counts)
  for (p in c(snakemake@output$pca, snakemake@output$distance,
              snakemake@output$volcano, snakemake@output$heatmap)) {
    png(p, width = 800, height = 600, res = 150)
    plot.new(); title("No viral genes detected in count matrix")
    dev.off()
  }
  message("Placeholder outputs created. Pipeline can continue.")
  teardown_log(log_con)
  quit(save = "no", status = 0)
}

# 3. Load data
dat <- load_common_samples(snakemake@input$counts, snakemake@input$samples)
if (nrow(dat$counts) < 2) {
  warning("Very few viral genes detected. DESeq2 analysis may be limited.")
}

# 4. Design and config
cfg <- setup_design(dat$metadata, snakemake@config)

# 5. Pre-filter (more permissive defaults for sparse viral counts)
counts_filt <- prefilter_counts(dat$counts, snakemake@config,
                                default_min_count = 5, default_min_samples = 2)
if (nrow(counts_filt) == 0) {
  message("All viral genes filtered out. Creating placeholder outputs.")
  empty_df <- data.frame(
    GeneID = character(), GeneSymbol = character(), GeneBiotype = character(),
    baseMean = numeric(), log2FoldChange = numeric(), lfcSE = numeric(),
    stat = numeric(), pvalue = numeric(), padj = numeric())
  write.table(empty_df, snakemake@output$results, sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(empty_df, snakemake@output$results_filtered, sep = "\t",
              quote = FALSE, row.names = FALSE)
  writeLines("gene_id", snakemake@output$normalized_counts)
  for (p in c(snakemake@output$pca, snakemake@output$distance,
              snakemake@output$volcano, snakemake@output$heatmap)) {
    png(p, width = 800, height = 600, res = 150)
    plot.new(); title("No viral genes passed pre-filtering")
    dev.off()
  }
  message("Placeholder outputs created. Pipeline can continue.")
  teardown_log(log_con)
  quit(save = "no", status = 0)
}

# 6. Build DESeq2 dataset manually (poscounts size factors needed)
dds <- DESeqDataSetFromMatrix(countData = round(counts_filt),
                               colData   = cfg$metadata,
                               design    = cfg$design_formula)
cnts <- counts(dds)
sf <- DESeq2::estimateSizeFactorsForMatrix(cnts, type = "poscounts")
if (any(is.na(sf))) {
  message("Some size factors were NA (e.g. zero viral counts); setting to 1.")
  sf[is.na(sf)] <- 1
}
sizeFactors(dds) <- sf
dds <- DESeq(dds)

# 7. Extract results
res <- extract_results(dds, cfg$metadata, snakemake@config)
res_df <- as.data.frame(res)
res_df$GeneID     <- rownames(res_df)
res_df$GeneSymbol <- res_df$GeneID
res_df$GeneBiotype <- "viral_gene"

# 8. Write outputs
write_outputs(res_df, dds, cfg$padj_threshold, snakemake@output)

# 9. Transform for visualization (VST with log2 fallback)
n_genes <- nrow(dds)
vsd <- tryCatch(vst(dds, blind = FALSE, nsub = min(1000L, n_genes)),
                error = function(e) NULL)
use_vst <- !is.null(vsd)
if (!use_vst) {
  message("VST failed (typical for few viral genes); using log2(norm + 1).")
  norm_counts <- counts(dds, normalized = TRUE)
  transform_data <- log2(norm_counts + 1)
} else {
  transform_data <- vsd
}

# 10. Visualizations
make_visualizations(transform_data, res_df, res, cfg$metadata,
                     cfg$condition_col, cfg$padj_threshold,
                     gene_label_col = "GeneID",
                     output = snakemake@output,
                     title_text = " (Viral Genes)",
                     use_vst = use_vst,
                     label_threshold = cfg$padj_threshold,
                     volcano_alpha = 0.6)

message("Viral gene DESeq2 analysis completed successfully")

# 11. Teardown
teardown_log(log_con)
