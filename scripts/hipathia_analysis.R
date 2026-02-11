# Redirect output to a log file if Snakemake is capturing it
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Install HiPathia from Bioconductor if not available (conda version may be missing pathway-viewer)
if (!require("hipathia", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("hipathia", ask = FALSE, update = FALSE)
    library(hipathia)
} else {
    # Check if pathway-viewer exists - if not, reinstall from Bioconductor
    hipathia_home_check <- system.file(package = "hipathia")
    pathway_viewer_check <- file.path(hipathia_home_check, "pathway-viewer")
    if (!dir.exists(pathway_viewer_check)) {
        message("HiPathia found but pathway-viewer is missing. Reinstalling from Bioconductor...")
        if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("hipathia", ask = FALSE, update = TRUE)
        # Reload to get the new installation
        detach("package:hipathia", unload = TRUE)
        library(hipathia)
    }
}

library(ggplot2)
library(pheatmap)

message("Loading normalized counts and metadata...")
# Load normalized counts
normalized_counts <- read.table(snakemake@input$normalized_counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.table(snakemake@input$samples, header=TRUE, row.names=1, sep="\t")

# Find common samples
common_samples <- intersect(colnames(normalized_counts), rownames(metadata))
if (length(common_samples) == 0) {
    stop("No common samples found between normalized counts and metadata file!")
}

# Filter to common samples
normalized_counts <- normalized_counts[, common_samples, drop=FALSE]
metadata <- metadata[common_samples, , drop=FALSE]

message(paste("Analyzing", ncol(normalized_counts), "samples with", nrow(normalized_counts), "genes"))

# Remove version numbers from Ensembl IDs for HiPathia compatibility
rownames(normalized_counts) <- sub("\\.[0-9]+$", "", rownames(normalized_counts))

message("Loading pathway data from HiPathia...")
# Load pathway data (this will download if not available)
pathways <- load_pathways("hsa")  # hsa = Homo sapiens

message("Preparing expression matrix for HiPathia...")
# Convert to matrix (HiPathia expects a matrix)
expr_matrix <- as.matrix(normalized_counts)

message("Computing pathway activity...")
# Compute pathway activity (returns subpathway-level activity)
# This is HiPathia's core function that computes signal transduction
pathway_activity <- hipathia(expr_matrix, pathways, decompose = FALSE, verbose = TRUE)

# Extract subpathway activity values using HiPathia's function
subpathway_values <- get_paths_data(pathway_activity)
message(paste("Extracted", nrow(subpathway_values), "subpathways"))
message(paste("Sample subpathway names:", paste(head(rownames(subpathway_values), 5), collapse=", ")))

# Normalize pathway values by rows (HiPathia's normalization function)
# normalize_paths may return a SummarizedExperiment, so extract the assay data
subpathway_normalized_obj <- normalize_paths(subpathway_values, pathways)
# Check if it's a SummarizedExperiment or matrix
if (class(subpathway_normalized_obj)[1] == "SummarizedExperiment") {
    subpathway_values_normalized <- assay(subpathway_normalized_obj)
} else {
    subpathway_values_normalized <- as.matrix(subpathway_normalized_obj)
}

# Get pathway list
pathway_list <- get_pathways_list(pathways)
message(paste("Found", length(pathway_list), "pathways"))

# Aggregate subpathways to pathway level
# Extract pathway ID from subpathway names
# Subpathway names are in format "P-hsa04010-1" or "P-hsa04010-2" etc.
# Extract the pathway ID (hsa04010) from this format
subpathway_to_pathway <- sub("^P-(hsa[0-9]+)-.*$", "\\1", rownames(subpathway_values_normalized))
# If extraction failed (doesn't match pattern), try alternative patterns
if (any(subpathway_to_pathway == rownames(subpathway_values_normalized))) {
    # Try other patterns
    subpathway_to_pathway <- sub("^.*(hsa[0-9]+).*$", "\\1", rownames(subpathway_values_normalized))
}

pathway_values_list <- list()
for (pathway_id in pathway_list) {
    # Find all subpathways belonging to this pathway
    matching_indices <- which(subpathway_to_pathway == pathway_id)
    if (length(matching_indices) > 0) {
        matching_subpaths <- rownames(subpathway_values_normalized)[matching_indices]
        # Calculate mean activity across subpathways for this pathway
        pathway_values_list[[pathway_id]] <- colMeans(subpathway_values_normalized[matching_subpaths, , drop=FALSE])
    }
}

# Convert to matrix
if (length(pathway_values_list) > 0) {
    pathway_values <- do.call(rbind, pathway_values_list)
    rownames(pathway_values) <- names(pathway_values_list)
    message(paste("Successfully aggregated", length(pathway_values_list), "pathways from subpathways"))
} else {
    # Fallback: if aggregation fails, use subpathways directly but extract pathway IDs
    warning("Could not aggregate to pathway level. Using subpathway data with pathway IDs.")
    pathway_values <- subpathway_values_normalized
    # Extract pathway ID from subpathway names for labeling
    pathway_ids <- sub("[:_].*$", "", rownames(pathway_values))
    message(paste("Using subpathway-level data. Pathway IDs extracted from", length(unique(pathway_ids)), "unique pathways"))
}

# Convert to data frame for easier handling
pathway_df <- as.data.frame(pathway_values)
pathway_df$Pathway <- rownames(pathway_df)

# Get pathway names for better readability (with error handling)
# Try to get names one by one to handle unrecognized pathway IDs gracefully
pathway_names <- sapply(rownames(pathway_df), function(pw_id) {
    tryCatch({
        name <- get_path_names(pathways, pw_id)
        if (length(name) > 0 && !is.na(name) && name != "") {
            name
        } else {
            pw_id  # Return ID if name is empty
        }
    }, error = function(e) {
        pw_id  # Return ID if name lookup fails
    })
})

pathway_df$PathwayName <- pathway_names

# Reorder columns
pathway_df <- pathway_df[, c("Pathway", "PathwayName", setdiff(colnames(pathway_df), c("Pathway", "PathwayName")))]

# Write pathway activity results
write.table(pathway_df, snakemake@output$pathway_activity, sep="\t", quote=FALSE, row.names=FALSE)
message(paste("Saved pathway activity for", nrow(pathway_df), "pathways"))

# Create pathway comparison plots using HiPathia's visualization
message("Creating pathway comparison plots...")

# Get condition/group labels for comparison
if ("Condition" %in% colnames(metadata)) {
    sample_groups <- as.character(metadata$Condition)
    group_col <- "Condition"
} else if ("Control" %in% colnames(metadata)) {
    sample_groups <- as.character(metadata$Control)
    group_col <- "Control"
} else {
    warning("No Condition or Control column found. Cannot create comparison plots.")
    sample_groups <- NULL
}

if (!is.null(sample_groups) && length(unique(sample_groups)) >= 2) {
    # Perform Wilcoxon test to compare pathways between groups
    # do_wilcoxon requires g1 and g2 parameters to specify which groups to compare
    message("Performing Wilcoxon test for pathway comparisons...")
    
    # Get unique groups
    unique_groups <- unique(sample_groups)
    message(paste("Found groups:", paste(unique_groups, collapse=", ")))
    
    # Use the normalized subpathway values matrix directly
    # Ensure sample groups are in the same order as matrix columns
    sample_groups_ordered <- sample_groups[match(colnames(subpathway_values_normalized), rownames(metadata))]
    
    # Determine which groups to compare (typically first two groups, or Experiment vs Negative/Control)
    if ("Experiment" %in% unique_groups && "Negative" %in% unique_groups) {
        g1 <- "Experiment"
        g2 <- "Negative"
    } else if ("Experiment" %in% unique_groups && "Control" %in% unique_groups) {
        g1 <- "Experiment"
        g2 <- "Control"
    } else {
        # Use first two groups found
        g1 <- unique_groups[1]
        g2 <- unique_groups[2]
    }
    
    message(paste("Comparing", g1, "vs", g2))
    
    # Perform Wilcoxon test using the matrix with specified groups
    comp <- do_wilcoxon(subpathway_values_normalized, sample_groups_ordered, g1 = g1, g2 = g2)
    
    # Create heatmap using HiPathia's built-in heatmap_plot function
    message("Creating subpathway activity heatmap...")
    message(paste("Matrix dimensions:", nrow(subpathway_values_normalized), "subpathways x", ncol(subpathway_values_normalized), "samples"))
    
    # Check for problematic values (NA, NaN, Inf)
    has_na <- any(is.na(subpathway_values_normalized))
    has_nan <- any(is.nan(subpathway_values_normalized))
    has_inf <- any(is.infinite(subpathway_values_normalized))
    
    if (has_na || has_nan || has_inf) {
        message("Found problematic values. Cleaning data...")
        # Replace problematic values with 0 or row means
        subpathway_values_clean <- subpathway_values_normalized
        subpathway_values_clean[is.na(subpathway_values_clean) | is.nan(subpathway_values_clean) | is.infinite(subpathway_values_clean)] <- 0
        subpathway_values_normalized <- subpathway_values_clean
        message("Data cleaned")
    }
    
    # Use HiPathia's heatmap_plot function
    # This creates a proper heatmap showing subpathway activity across samples
    # Start with minimal clustering to avoid errors
    png(snakemake@output$pathway_heatmap, width=2000, height=min(3000, 200 + nrow(subpathway_values_normalized) * 2), res=150)
    
    tryCatch({
        # Try with sample clustering only (no variable clustering to avoid removing rows)
        heatmap_plot(subpathway_values_normalized, 
                     group = sample_groups_ordered,
                     colors = "hipathia",  # HiPathia color scheme: Green (low) -> White -> Orange (high)
                     sample_clust = TRUE,   # Cluster samples
                     variable_clust = FALSE, # Disable variable clustering to avoid row removal
                     scale = TRUE,          # Scale rows to [0,1]
                     main = "HiPathia Subpathway Activity Heatmap")
        message("Subpathway activity heatmap saved")
    }, error = function(e) {
        warning(paste("Error creating heatmap with clustering:", e$message))
        # Try without any clustering
        tryCatch({
            heatmap_plot(subpathway_values_normalized, 
                         group = sample_groups_ordered,
                         colors = "hipathia",
                         sample_clust = FALSE,  # Disable sample clustering
                         variable_clust = FALSE, # Disable variable clustering
                         scale = TRUE,
                         main = "HiPathia Subpathway Activity Heatmap")
            message("Subpathway activity heatmap saved (without clustering)")
        }, error = function(e2) {
            warning(paste("Error creating heatmap even without clustering:", e2$message))
            # Last resort: simple plot
            plot.new()
            text(0.5, 0.5, paste("Could not create heatmap.\nError:", e2$message), cex = 1.2)
        })
    })
    
    dev.off()
} else {
    warning("Insufficient groups for comparison. Creating placeholder plot.")
    png(snakemake@output$pathway_heatmap, width=800, height=600, res=150)
    plot.new()
    text(0.5, 0.5, "Pathway comparison requires at least 2 groups", cex = 1.5)
    dev.off()
}

# Create interactive HTML report using HiPathia's built-in functions
message("Creating interactive HiPathia report...")

# Use the same sample_groups from above (already defined)
if (!is.null(sample_groups) && length(unique(sample_groups)) >= 2) {
    # Perform Wilcoxon test if not already done
    if (!exists("comp")) {
        message("Performing Wilcoxon test for pathway comparisons...")
        # Get unique groups and determine comparison groups
        unique_groups <- unique(sample_groups)
        if ("Experiment" %in% unique_groups && "Negative" %in% unique_groups) {
            g1 <- "Experiment"
            g2 <- "Negative"
        } else if ("Experiment" %in% unique_groups && "Control" %in% unique_groups) {
            g1 <- "Experiment"
            g2 <- "Control"
        } else {
            g1 <- unique_groups[1]
            g2 <- unique_groups[2]
        }
        sample_groups_ordered <- sample_groups[match(colnames(subpathway_values_normalized), rownames(metadata))]
        comp <- do_wilcoxon(subpathway_values_normalized, sample_groups_ordered, g1 = g1, g2 = g2)
    }
    
    # Create interactive HTML report
    # Use relative path (simpler and more reliable)
    report_folder <- dirname(snakemake@output$report_html)
    
    message(paste("Creating report in:", report_folder))
    message(paste("Current working directory:", getwd()))
    
    # Create directory using relative path (should work from project root)
    if (!dir.exists(report_folder)) {
        message(paste("Creating report directory:", report_folder))
        success <- dir.create(report_folder, recursive = TRUE, showWarnings = TRUE)
        if (!success) {
            warning("dir.create() may have failed, but continuing...")
        }
    }
    
    # Verify directory exists (try both relative and check if it was created)
    if (!dir.exists(report_folder)) {
        # Try creating parent directory first
        parent <- dirname(report_folder)
        if (!dir.exists(parent)) {
            dir.create(parent, recursive = TRUE, showWarnings = TRUE)
        }
        dir.create(report_folder, recursive = TRUE, showWarnings = TRUE)
    }
    
    # Use HiPathia's create_report function
    message("Calling HiPathia's create_report()...")
    
    # Check if HiPathia's pathway-viewer files are available
    # create_report needs these files to create the interactive report
    hipathia_home <- system.file(package = "hipathia")
    pathway_viewer_path <- file.path(hipathia_home, "pathway-viewer")
    message(paste("HiPathia package location:", hipathia_home))
    message(paste("Pathway-viewer path:", pathway_viewer_path))
    message(paste("Pathway-viewer exists:", dir.exists(pathway_viewer_path)))
    
    # Check for pathway-viewer in various possible locations
    alt_paths <- c(
        pathway_viewer_path,
        file.path(hipathia_home, "extdata", "pathway-viewer"),
        file.path(hipathia_home, "inst", "pathway-viewer"),
        file.path(hipathia_home, "..", "hipathia", "pathway-viewer")
    )
    
    pathway_viewer_found <- FALSE
    for (alt_path in alt_paths) {
        if (dir.exists(alt_path)) {
            message(paste("Found pathway-viewer at:", alt_path))
            pathway_viewer_path <- alt_path
            pathway_viewer_found <- TRUE
            break
        }
    }
    
    if (!pathway_viewer_found) {
        warning("WARNING: pathway-viewer directory not found in HiPathia installation!")
        warning("This is required for create_report() to work.")
        warning("\nThe conda installation of HiPathia may not include pathway-viewer files.")
        warning("create_report() will likely fail. See error message for solution.")
    }
    
    # Clean up any existing directory
    if (dir.exists(report_folder)) {
        message("Removing existing report directory...")
        unlink(report_folder, recursive = TRUE)
    }
    
    # Create the directory first (empty) - this might help with file.copy issues
    dir.create(report_folder, recursive = TRUE, showWarnings = FALSE)
    
    message(paste("Creating report in:", report_folder))
    message(paste("Current working directory:", getwd()))
    
    # Try creating the report
    # The directory exists but is empty, which should help with file.copy
    tryCatch({
        report_result <- create_report(comp, pathways, report_folder)
        message("HiPathia report created successfully!")
        message(paste("Report location:", report_result))
        
        # Verify index.html was created
        index_file <- file.path(report_result, "index.html")
        if (file.exists(index_file)) {
            message(paste("Report index file created:", index_file))
        } else {
            # Check what was actually created
            if (dir.exists(report_result)) {
                files_created <- list.files(report_result, recursive = TRUE)
                message(paste("Files created in report folder:", paste(files_created, collapse=", ")))
            }
        }
        
        message("\nTo view the report:")
        message(paste("  1. Open", index_file, "in a web browser"))
        message(paste("  2. Or use in R: visualize_report('", report_result, "')", sep=""))
        
    }, error = function(e) {
        error_msg <- e$message
        warning(paste("Error creating HiPathia report:", error_msg))
        
        # Check if it's the pathway-viewer issue
        if (grepl("more 'from' files than 'to' files", error_msg) || !dir.exists(pathway_viewer_path)) {
            message("\n" , paste(rep("=", 60), collapse=""))
            message("DIAGNOSIS: pathway-viewer component is missing from HiPathia installation")
            message(paste(rep("=", 60), collapse=""))
            message("\nThe conda-installed HiPathia package is missing the 'pathway-viewer' directory")
            message("which is required for create_report() to generate interactive reports.")
            message("\nSOLUTION:")
            message("1. Install HiPathia from Bioconductor instead of conda:")
            message("   In R, run:")
            message("     if (!require('BiocManager')) install.packages('BiocManager')")
            message("     BiocManager::install('hipathia', force = TRUE)")
            message("\n2. Then re-run this analysis")
            message("\n3. Alternatively, you can use the pathway activity TSV file and heatmap")
            message("   which are already generated successfully.")
            message(paste(rep("=", 60), collapse=""), "\n")
        }
        
        # Create a helpful HTML file explaining the issue
        index_content <- c(
            "<!DOCTYPE html>",
            "<html><head><title>HiPathia Report</title>",
            "<style>body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }",
            "pre { background: #f4f4f4; padding: 15px; border-radius: 5px; }",
            ".error { color: #d32f2f; background: #ffebee; padding: 15px; border-radius: 5px; }",
            ".solution { color: #1976d2; background: #e3f2fd; padding: 15px; border-radius: 5px; margin-top: 20px; }</style>",
            "</head><body>",
            "<h1>HiPathia Pathway Analysis Report</h1>",
            "<div class='error'>",
            "<h2>Interactive Report Generation Failed</h2>",
            paste("<p><strong>Error:</strong>", error_msg, "</p>"),
            "</div>",
            "<div class='solution'>",
            "<h2>Solution</h2>",
            "<p>The conda-installed HiPathia package is missing the 'pathway-viewer' component.</p>",
            "<p>To fix this, install HiPathia from Bioconductor in R:</p>",
            "<pre>",
            "if (!require('BiocManager')) install.packages('BiocManager')",
            "BiocManager::install('hipathia', force = TRUE)",
            "</pre>",
            "<p>Then re-run the analysis.</p>",
            "</div>",
            "<h2>Available Results</h2>",
            "<ul>",
            "<li><strong>Pathway Activity Data:</strong> <code>results/hipathia_pathway_activity.tsv</code></li>",
            "<li><strong>Pathway Heatmap:</strong> <code>results/hipathia_pathway_heatmap.png</code></li>",
            "</ul>",
            "</body></html>"
        )
        writeLines(index_content, snakemake@output$report_html)
        message("Created HTML file with diagnostic information and solution.")
    })
    
    message(paste("Interactive report created in:", report_folder))
    message("To view the report, open the index.html file in a web browser")
    message(paste("Or use: visualize_report('", report_folder, "') in R", sep=""))
} else {
    warning("Insufficient groups for comparison. Skipping interactive report creation.")
    # Create empty directory and placeholder file for output
    report_folder <- normalizePath(dirname(snakemake@output$report_html), mustWork = FALSE)
    dir.create(report_folder, recursive = TRUE, showWarnings = FALSE)
    writeLines("Interactive report requires at least 2 groups for comparison.", 
               file.path(report_folder, "README.txt"))
    # Create placeholder index.html
    writeLines("<html><body><h1>HiPathia Report</h1><p>Interactive report requires at least 2 groups for comparison.</p></body></html>",
               snakemake@output$report_html)
}

message("HiPathia analysis completed successfully!")
