# DESeq2 Analysis — combined / human counts
source("scripts/deseq2_core.R")
library(biomaRt)

# 1. Setup
log_con <- setup_log(snakemake@log[[1]])

# 2. Load data
dat <- load_common_samples(snakemake@input$counts, snakemake@input$samples)

# 3. Design and config
cfg <- setup_design(dat$metadata, snakemake@config)

# 4. Pre-filter
counts_filt <- prefilter_counts(dat$counts, snakemake@config,
                                default_min_count = 10, default_min_samples = 3)

# 5. Run DESeq2
dds <- run_deseq2(counts_filt, cfg$metadata, cfg$design_formula)

# 6. Extract results
res <- extract_results(dds, cfg$metadata, snakemake@config)
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

# 7a. Map Ensembl IDs to gene symbols (human-specific)
message("Mapping Ensembl IDs to gene symbols...")
res_df$GeneSymbol <- NA
res_df$GeneBiotype <- NA
ensembl_ids <- sub("\\.[0-9]+$", "", res_df$GeneID)
tryCatch({
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  batch_size <- 10000
  gene_map_list <- list()
  for (i in seq(1, length(ensembl_ids), by = batch_size)) {
    end_idx <- min(i + batch_size - 1, length(ensembl_ids))
    gene_map_batch <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
      filters     = "ensembl_gene_id",
      values      = ensembl_ids[i:end_idx],
      mart        = ensembl)
    gene_map_list[[length(gene_map_list) + 1]] <- gene_map_batch
  }
  gene_map <- do.call(rbind, gene_map_list)
  gene_map_unique <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  res_df$GeneSymbol <- gene_map_unique$hgnc_symbol[
    match(ensembl_ids, gene_map_unique$ensembl_gene_id)]
  res_df$GeneBiotype <- gene_map_unique$gene_biotype[
    match(ensembl_ids, gene_map_unique$ensembl_gene_id)]
  res_df$GeneSymbol[res_df$GeneSymbol == ""] <- NA
  message(paste("Successfully mapped", sum(!is.na(res_df$GeneSymbol)),
                "genes to symbols"))
}, error = function(e) {
  warning(paste("Failed to map gene symbols via biomaRt:", e$message))
  warning("Results will contain only Ensembl IDs")
})

# 7b. Write outputs
write_outputs(res_df, dds, cfg$padj_threshold, snakemake@output)

# 8. VST and visualizations
vsd <- vst(dds, blind = FALSE)
make_visualizations(vsd, res_df, res, cfg$metadata, cfg$condition_col,
                     cfg$padj_threshold, gene_label_col = "GeneSymbol",
                     output = snakemake@output, title_text = "",
                     use_vst = TRUE, label_threshold = 0.000001)

# 9. Teardown
teardown_log(log_con)
