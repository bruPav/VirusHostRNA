# =============================================================================
# Generate a self-contained HTML summary report for the VirusHostRNA pipeline.
# Uses only base R — no RMarkdown, no pandoc, no extra packages.
# =============================================================================
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

# --- Read DE gene counts from human results table ----------------------------
de_human <- read.table(snakemake@input$deg_human, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
n_de_human <- sum(!is.na(de_human$padj) & de_human$padj < 0.05, na.rm = TRUE)
n_total_human <- nrow(de_human)

# --- Read DE gene counts from viral results table (may be empty) -------------
de_viral <- read.table(snakemake@input$deg_viral, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
n_de_viral <- sum(!is.na(de_viral$padj) & de_viral$padj < 0.05, na.rm = TRUE)
n_total_viral <- nrow(de_viral)

# --- Helper: relative link to a file in results/ -----------------------------
rel <- function(path) basename(path)

# --- Helper: embed an image as base64 data URI or use relative path ----------
img_tag <- function(src, alt = "", width = "100%") {
  sprintf('<img src="%s" alt="%s" style="max-width:%s; height:auto; display:block; margin:10px 0;"/>',
          src, alt, width)
}

# --- Helper: write a section block -------------------------------------------
section <- function(heading, ...) {
  paste0("<h2>", heading, "</h2>\n", paste(..., collapse = "\n"), "\n<hr/>\n")
}

# --- Helper: clickable file link ---------------------------------------------
file_link <- function(path, label) {
  f <- basename(path)
  if (file.exists(path)) {
    sprintf('<a href="%s">%s</a>', f, label)
  } else {
    sprintf('<span style="color:#999">[missing] %s</span>', label)
  }
}

# --- Build HTML --------------------------------------------------------------
html <- c(
  "<!DOCTYPE html>",
  "<html><head><meta charset='utf-8'/>",
  "<title>VirusHostRNA Pipeline Report</title>",
  "<style>",
  "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;",
  "       max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333;",
  "       line-height: 1.6; }",
  "h1 { border-bottom: 2px solid #2563eb; padding-bottom: 8px; color: #1e40af; }",
  "h2 { color: #2563eb; margin-top: 30px; }",
  ".card { border: 1px solid #e5e7eb; border-radius: 8px; padding: 16px;",
  "         margin: 12px 0; background: #f9fafb; }",
  ".stat { display: inline-block; min-width: 140px; margin: 4px 8px; }",
  ".stat-val { font-size: 1.3em; font-weight: bold; color: #2563eb; }",
  ".stat-label { font-size: 0.85em; color: #6b7280; }",
  "a { color: #2563eb; }",
  "table { border-collapse: collapse; width: 100%; }",
  "th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid #e5e7eb; }",
  "th { background: #f3f4f6; }",
  "</style></head><body>",
  "<h1>VirusHostRNA Pipeline Report</h1>",
  paste("<p><strong>Date:</strong>", Sys.Date(), "</p>"),
  paste("<p><strong>R version:</strong>", R.version.string, "</p>"),

  # ---- Quick Stats ----------------------------------------------------------
  section("Summary",
    "<div class='card'>",
    sprintf("<div class='stat'><div class='stat-val'>%d</div><div class='stat-label'>Human genes tested</div></div>", n_total_human),
    sprintf("<div class='stat'><div class='stat-val'>%d</div><div class='stat-label'>Human DE genes</div></div>", n_de_human),
    sprintf("<div class='stat'><div class='stat-val'>%d</div><div class='stat-label'>Viral genes tested</div></div>", n_total_viral),
    sprintf("<div class='stat'><div class='stat-val'>%d</div><div class='stat-label'>Viral DE genes</div></div>", n_de_viral),
    "</div>"
  ),

  # ---- QC Report ------------------------------------------------------------
  section("Quality Control",
    "<ul>",
    sprintf("<li>%s</li>", file_link(snakemake@input$multiqc, "MultiQC Report")),
    "</ul>"
  ),

  # ---- Human DESeq2 ---------------------------------------------------------
  section("Human Gene Expression",
    img_tag(rel(snakemake@input$pca_human), "Human PCA"),
    img_tag(rel(snakemake@input$volcano_human), "Human Volcano"),
    img_tag(rel(snakemake@input$heatmap_human), "Human Heatmap"),
    img_tag(rel(snakemake@input$distance_human), "Human Distance Matrix"),
    "<ul>",
    sprintf("<li>%s (%d genes)</li>",
            file_link(snakemake@input$deg_human, "Full DE results"),
            n_total_human),
    sprintf("<li>%s (%d significant)</li>",
            file_link(snakemake@input$deg_sig_human, "Significant DE genes"),
            n_de_human),
    sprintf("<li>%s</li>",
            file_link(snakemake@input$norm_human, "Normalized counts (human)")),
    "</ul>"
  ),

  # ---- Viral DESeq2 ---------------------------------------------------------
  section("Viral Gene Expression",
    img_tag(rel(snakemake@input$pca_viral), "Viral PCA"),
    img_tag(rel(snakemake@input$volcano_viral), "Viral Volcano"),
    img_tag(rel(snakemake@input$heatmap_viral), "Viral Heatmap"),
    img_tag(rel(snakemake@input$distance_viral), "Viral Distance Matrix"),
    "<ul>",
    sprintf("<li>%s (%d genes)</li>",
            file_link(snakemake@input$deg_viral, "Full DE results"),
            n_total_viral),
    sprintf("<li>%s (%d significant)</li>",
            file_link(snakemake@input$deg_sig_viral, "Significant DE genes"),
            n_de_viral),
    sprintf("<li>%s</li>",
            file_link(snakemake@input$norm_viral, "Normalized counts (viral)")),
    "</ul>"
  ),

  # ---- Enrichment -----------------------------------------------------------
  section("Functional Enrichment (Human)",
    img_tag(rel(snakemake@input$enrich_kegg_dot), "KEGG GSEA"),
    img_tag(rel(snakemake@input$enrich_gobp_dot), "GO BP GSEA"),
    img_tag(rel(snakemake@input$enrich_kegg_ora_dot), "KEGG ORA"),
    img_tag(rel(snakemake@input$enrich_gobp_ora_dot), "GO BP ORA"),
    "<ul>",
    sprintf("<li>%s</li>", file_link(snakemake@input$enrich_kegg, "KEGG GSEA table")),
    sprintf("<li>%s</li>", file_link(snakemake@input$enrich_gobp, "GO BP GSEA table")),
    sprintf("<li>%s</li>", file_link(snakemake@input$enrich_kegg_ora, "KEGG ORA table")),
    sprintf("<li>%s</li>", file_link(snakemake@input$enrich_gobp_ora, "GO BP ORA table")),
    "</ul>"
  ),

  # ---- All Output Files -----------------------------------------------------
  section("All Output Files",
    "<table><tr><th>File</th><th>Description</th></tr>",
    sprintf("<tr><td>%s</td><td>Count matrix (combined human + viral)</td></tr>",
            file_link(snakemake@input$count_matrix, "gene_counts_matrix.tsv")),
    sprintf("<tr><td>%s</td><td>Human count matrix</td></tr>",
            file_link(snakemake@input$count_human, "gene_counts_matrix_human.tsv")),
    sprintf("<tr><td>%s</td><td>Viral count matrix</td></tr>",
            file_link(snakemake@input$count_viral, "gene_counts_matrix_viral.tsv")),
    sprintf("<tr><td>%s</td><td>Combined DESeq2 results</td></tr>",
            file_link(snakemake@input$deg_combined, "deg_results.tsv")),
    sprintf("<tr><td>%s</td><td>Normalized counts (combined)</td></tr>",
            file_link(snakemake@input$norm_combined, "normalized_counts.tsv")),
    sprintf("<tr><td>%s</td><td>Strandedness summary</td></tr>",
            file_link(snakemake@input$strandedness, "strandedness_summary.tsv")),
    sprintf("<tr><td>%s</td><td>MultiQC report</td></tr>",
            file_link(snakemake@input$multiqc, "multiqc_report.html")),
    "</table>"
  ),

  "<hr/>",
  "<footer style='color:#9ca3af; font-size:0.85em; margin-top:40px;'>",
  "Generated by VirusHostRNA pipeline",
  paste(" —", Sys.time()),
  "</footer>",
  "</body></html>"
)

# --- Write output ------------------------------------------------------------
writeLines(html, snakemake@output$report)
message(paste("HTML report written to", snakemake@output$report))

sink(type = "message")
sink()
close(log)
