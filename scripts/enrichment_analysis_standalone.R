#!/usr/bin/env Rscript
# Pathway and GO enrichment with clusterProfiler (ORA + GSEA)
# Uses DESeq2 results (human genes). Run from project root.
#
# Use the dedicated enrichment env (pre-built clusterProfiler from conda):
#   conda env create -f envs/enrichment.yaml   # once
#   conda activate enrichment && Rscript scripts/enrichment_analysis_standalone.R

# --- Paths ---
# Detect project root: prefer working directory, then look for known marker files
# in parent directories. No hardcoded absolute paths.
find_project_root <- function() {
  # If running from project root, results/ should be present or creatable
  if (file.exists("Snakefile") || file.exists("config.yaml") || file.exists("results/deg_results_human.tsv")) {
    return(getwd())
  }
  # Walk up directory tree looking for Snakefile
  d <- getwd()
  for (i in 1:5) {
    d <- dirname(d)
    if (file.exists(file.path(d, "Snakefile"))) return(d)
  }
  return(getwd())  # fallback
}
project_dir <- find_project_root()

deseq2_file <- file.path(project_dir, "results", "deg_results_human.tsv")
out_dir     <- file.path(project_dir, "results", "enrichment")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Enrichment analysis (ORA + GSEA). First run may take 1–2 min (KEGG/GO API).")

# Load packages (install if missing; enrichplot optional - script works without it)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
if (has_enrichplot) library(enrichplot)

message("Loading DESeq2 results...")
res <- read.table(deseq2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Strip Ensembl version for mapping
res$Ensembl <- sub("\\.[0-9]+$", "", res$GeneID)

# Map Ensembl -> Entrez for KEGG/GO (org.Hs.eg.db)
# Many IDs "fail to map" by design: lncRNAs, pseudogenes, and some genes are not in the DB.
# KEGG/GO enrichment uses the successfully mapped set; this is standard.
message("Mapping Ensembl to Entrez IDs (unmapped: often lncRNA/pseudogenes or not in DB)...")
id_map <- suppressWarnings(bitr(unique(res$Ensembl), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
# If 1:many Ensembl->Entrez, keep first Entrez per Ensembl
id_map <- id_map[!duplicated(id_map$ENSEMBL), ]
res$Entrez <- id_map$ENTREZID[match(res$Ensembl, id_map$ENSEMBL)]
n_mapped <- sum(!is.na(res$Entrez))
n_total  <- length(unique(res$Ensembl))
message(paste("  Mapped", n_mapped, "of", n_total, "unique Ensembl IDs to Entrez (", round(100 * n_mapped / n_total, 1), "%).", sep = ""))
res <- res[!is.na(res$Entrez) & !duplicated(res$Entrez), ]

# Remove genes with NA stats (for ranking)
res_rank <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]
message(paste("Using", nrow(res_rank), "genes with valid stats for enrichment"))

# --- ORA: Over-representation analysis ---
# DE genes at padj < 0.05 (relax to 0.1 if very few genes)
padj_cutoff <- 0.05
de_entrez <- res_rank$Entrez[!is.na(res_rank$padj) & res_rank$padj < padj_cutoff]
if (length(de_entrez) < 10) {
  padj_cutoff <- 0.1
  de_entrez   <- res_rank$Entrez[!is.na(res_rank$padj) & res_rank$padj < padj_cutoff]
  message(paste("Few DE at 0.05; using padj <", padj_cutoff, "->", length(de_entrez), "genes"))
} else {
  message(paste("ORA: using", length(de_entrez), "DE genes (padj <", padj_cutoff, ")"))
}

ora_kegg <- NULL
ora_go   <- NULL
if (length(de_entrez) >= 5) {
  message("Running KEGG over-representation...")
  ora_kegg <- tryCatch({
    enrichKEGG(gene = de_entrez, organism = "hsa", pvalueCutoff = 0.2, qvalueCutoff = 0.2)
  }, error = function(e) { message("enrichKEGG: ", e$message); NULL })
  if (!is.null(ora_kegg) && nrow(ora_kegg@result) > 0) {
    write.table(ora_kegg@result, file.path(out_dir, "enrich_KEGG_ORA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("  KEGG ORA: ", nrow(ora_kegg@result), "terms -> enrich_KEGG_ORA.tsv"))
    tryCatch({
      if (has_enrichplot) {
        p <- dotplot(ora_kegg, showCategory = 15) + ggtitle("KEGG over-representation (DE genes)") +
          theme(axis.text.y = element_text(size = 7))
      } else {
        r <- head(ora_kegg@result[order(ora_kegg@result$pvalue), ], 15)
        r$Description <- substr(r$Description, 1, 50)
        p <- ggplot(r, aes(x = reorder(Description, -pvalue), y = -log10(pvalue), fill = p.adjust)) +
          geom_col() + coord_flip() + theme_minimal() + labs(x = "", title = "KEGG over-representation (DE genes)") +
          theme(axis.text.y = element_text(size = 7))
      }
      ggsave(file.path(out_dir, "enrich_KEGG_ORA_dotplot.png"), p, width = 8, height = max(5, min(15, nrow(ora_kegg@result)) * 0.3), dpi = 150)
    }, error = function(e) message("  Dotplot skipped: ", e$message))
  }
  message("Running GO BP over-representation...")
  ora_go <- tryCatch({
    enrichGO(gene = de_entrez, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2, qvalueCutoff = 0.2)
  }, error = function(e) { message("enrichGO: ", e$message); NULL })
  if (!is.null(ora_go) && nrow(ora_go@result) > 0) {
    write.table(ora_go@result, file.path(out_dir, "enrich_GOBP_ORA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("  GO BP ORA:", nrow(ora_go@result), "terms -> enrich_GOBP_ORA.tsv"))
    tryCatch({
      if (has_enrichplot) {
        p <- dotplot(ora_go, showCategory = 15) + ggtitle("GO BP over-representation (DE genes)") +
          theme(axis.text.y = element_text(size = 7))
      } else {
        r <- head(ora_go@result[order(ora_go@result$pvalue), ], 15)
        r$Description <- substr(r$Description, 1, 60)
        p <- ggplot(r, aes(x = reorder(Description, -pvalue), y = -log10(pvalue), fill = p.adjust)) +
          geom_col() + coord_flip() + theme_minimal() + labs(x = "", title = "GO BP over-representation (DE genes)") +
          theme(axis.text.y = element_text(size = 7))
      }
      ggsave(file.path(out_dir, "enrich_GOBP_ORA_dotplot.png"), p, width = 8, height = max(5, min(15, nrow(ora_go@result)) * 0.3), dpi = 150)
    }, error = function(e) message("  Dotplot skipped: ", e$message))
  }
} else {
  message("ORA skipped: fewer than 5 DE genes. Consider relaxing padj or running GSEA only.")
}

# --- GSEA: Gene Set Enrichment Analysis (ranked by signed statistic) ---
# Rank: sign(log2FC) * -log10(pvalue) so upregulated = positive, down = negative
res_rank$rank_stat <- sign(res_rank$log2FoldChange) * (-log10(res_rank$pvalue + 1e-300))
res_rank <- res_rank[is.finite(res_rank$rank_stat), ]
gene_list <- setNames(res_rank$rank_stat, res_rank$Entrez)
gene_list <- sort(gene_list, decreasing = TRUE)

message("Running KEGG GSEA...")
gsea_kegg <- tryCatch({
  gseKEGG(geneList = gene_list, organism = "hsa", pvalueCutoff = 0.2, pAdjustMethod = "BH")
}, error = function(e) { message("gseKEGG: ", e$message); NULL })
if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  write.table(gsea_kegg@result, file.path(out_dir, "enrich_KEGG_GSEA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste("  KEGG GSEA:", nrow(gsea_kegg@result), "terms -> enrich_KEGG_GSEA.tsv"))
  tryCatch({
    if (has_enrichplot) {
      p <- dotplot(gsea_kegg, showCategory = 15) + ggtitle("KEGG GSEA (ranked by signed -log10(p)") +
        theme(axis.text.y = element_text(size = 7))
    } else {
      r <- head(gsea_kegg@result[order(gsea_kegg@result$pvalue), ], 15)
      r$Description <- substr(r$Description, 1, 50)
      p <- ggplot(r, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
        geom_col() + coord_flip() + theme_minimal() + labs(x = "", title = "KEGG GSEA") +
        theme(axis.text.y = element_text(size = 7))
    }
    ggsave(file.path(out_dir, "enrich_KEGG_GSEA_dotplot.png"), p, width = 8, height = max(5, min(15, nrow(gsea_kegg@result)) * 0.25), dpi = 150)
  }, error = function(e) message("  Dotplot skipped: ", e$message))
} else {
  # Create empty output files so Snakemake doesn't fail on missing outputs
  message("  KEGG GSEA: no significant terms found. Creating empty output files.")
  writeLines("No KEGG GSEA results at the current p-value cutoff.", file.path(out_dir, "enrich_KEGG_GSEA.tsv"))
  png(file.path(out_dir, "enrich_KEGG_GSEA_dotplot.png"), width = 800, height = 600, res = 150)
  plot.new(); title("KEGG GSEA: no significant terms found")
  dev.off()
}

message("Running GO BP GSEA...")
gsea_go <- tryCatch({
  gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2, pAdjustMethod = "BH")
}, error = function(e) { message("gseGO: ", e$message); NULL })
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  write.table(gsea_go@result, file.path(out_dir, "enrich_GOBP_GSEA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste("  GO BP GSEA:", nrow(gsea_go@result), "terms -> enrich_GOBP_GSEA.tsv"))
  tryCatch({
    if (has_enrichplot) {
      p <- dotplot(gsea_go, showCategory = 15) + ggtitle("GO BP GSEA (ranked by signed -log10(p)") +
        theme(axis.text.y = element_text(size = 7))
    } else {
      r <- head(gsea_go@result[order(gsea_go@result$pvalue), ], 15)
      r$Description <- substr(r$Description, 1, 60)
      p <- ggplot(r, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
        geom_col() + coord_flip() + theme_minimal() + labs(x = "", title = "GO BP GSEA") +
        theme(axis.text.y = element_text(size = 7))
    }
    ggsave(file.path(out_dir, "enrich_GOBP_GSEA_dotplot.png"), p, width = 8, height = max(5, min(15, nrow(gsea_go@result)) * 0.3), dpi = 150)
  }, error = function(e) message("  Dotplot skipped: ", e$message))
} else {
  message("  GO BP GSEA: no significant terms found. Creating empty output files.")
  writeLines("No GO BP GSEA results at the current p-value cutoff.", file.path(out_dir, "enrich_GOBP_GSEA.tsv"))
  png(file.path(out_dir, "enrich_GOBP_GSEA_dotplot.png"), width = 800, height = 600, res = 150)
  plot.new(); title("GO BP GSEA: no significant terms found")
  dev.off()
}

message(paste("Enrichment results written to:", out_dir))
message("Done.")
