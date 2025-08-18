ogel$set('public', 'infercnv', function(
    gene_order_file,
    ref_group_names,
    celltype = NULL,
    cutoff = 0.1,
    out_dir = NULL,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = NULL,
    analysis_mode = "subclusters",
    assay_use = "RNA",
    data_use = "data"
) {
    
    self$say("Running InferCNV analysis")
    self$say(paste("Cutoff:", cutoff, "- genes with expression below this threshold in normal cells will be filtered"))
    self$say("Note: cutoff=0.1 is recommended for 10x data, cutoff=1 for Smart-seq2 data")
    
    # Check required parameters
    if (is.null(gene_order_file) || !file.exists(gene_order_file)) {
        stop("gene_order_file must be provided and must exist")
    }
    
    # Use provided celltype or default to self$celltype
    if (is.null(celltype)) {
        celltype <- self$celltype
    }
    
    if (is.null(celltype)) {
        stop("celltype must be specified either as parameter or in the ogel object")
    }
    
    if (is.null(ref_group_names)) {
        stop("ref_group_names must be specified")
    }
    
    # Ensure ref_group_names is character vector
    ref_group_names <- as.character(ref_group_names)
    
    # Remove any extra whitespace and ensure no special characters
    ref_group_names <- trimws(ref_group_names)
    
    # Get Seurat object
    obj <- self$get_data(data_use)
    
    if (!inherits(obj, 'Seurat')) {
        stop("InferCNV analysis requires a Seurat object")
    }
    
    # Set default parameters
    if (is.null(num_threads)) {
        num_threads <- self$threads
    }
    
    if (is.null(out_dir)) {
        out_dir <- file.path(self$path, "analysis", "infercnv")
    } else {
        if (dirname(out_dir) == ".") {
            # If it's just a folder name, put it in workspace/tag/analysis directory
            out_dir <- file.path(self$path, "analysis", out_dir)
        }
        # If it's a full path, use as is
    }
    
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
        self$say(paste("Created output directory:", out_dir))
    }
    
    # Load required library
    suppressMessages({
        if (!requireNamespace("infercnv", quietly = TRUE)) {
            stop("infercnv package is required but not installed")
        }
        library(infercnv)
    })
    
    verbose <- self$verbose
    
    # Prepare annotations file
    annotations_file <- file.path(out_dir, "cellAnnotations.txt")
    
    # Create cell annotations
    meta_data <- obj@meta.data
    annotations <- data.frame(
        cell_id = rownames(meta_data),
        group = trimws(as.character(meta_data[[celltype]])),  # Clean whitespace and ensure character
        stringsAsFactors = FALSE
    )
    
    # Get all unique groups from celltype
    all_groups <- unique(annotations$group)
    obs_group_names <- setdiff(all_groups, ref_group_names)
    
    # Ensure all group names are character vectors
    obs_group_names <- as.character(obs_group_names)
    
    self$say(paste("Using reference group:", paste(ref_group_names, collapse = ", ")))
    self$say(paste("Observation groups (auto-detected):", paste(obs_group_names, collapse = ", ")))
    
    # Log processed parameters for debugging
    self$say("Parameters being used:")
    self$say(paste("  gene_order_file:", gene_order_file))
    self$say(paste("  out_dir:", out_dir))
    self$say(paste("  celltype:", celltype))
    self$say(paste("  cutoff:", cutoff))
    self$say(paste("  cluster_by_groups:", cluster_by_groups))
    self$say(paste("  denoise:", denoise))
    self$say(paste("  HMM:", HMM))
    self$say(paste("  num_threads:", num_threads))
    self$say(paste("  analysis_mode:", analysis_mode))
    self$say(paste("  assay_use:", assay_use))
    
    # Check if ref_group_names exist in the data
    if (!all(ref_group_names %in% all_groups)) {
        missing_groups <- ref_group_names[!ref_group_names %in% all_groups]
        stop(paste("Reference groups not found in data:", paste(missing_groups, collapse = ", "), 
                  "\nAvailable groups:", paste(all_groups, collapse = ", ")))
    }
    
    write.table(annotations, file = annotations_file, 
               sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Get expression matrix
    expr_matrix <- GetAssayData(obj, assay = assay_use, layer = "counts")
    
    # Filter cells to match annotations
    common_cells <- intersect(colnames(expr_matrix), annotations$cell_id)
    expr_matrix <- expr_matrix[, common_cells]
    annotations <- annotations[annotations$cell_id %in% common_cells, ]
    
    # Store initial counts for summary
    initial_cells <- ncol(expr_matrix)
    initial_genes <- nrow(expr_matrix)
    
    self$say(paste("Input data:", initial_cells, "cells and", initial_genes, "genes"))
    
    # Create infercnv object
    self$say("Creating infercnv object...")
    
    # Try-catch to provide better error information
    tryCatch({
        # Suppress InferCNV verbose output
        suppressMessages({
            infercnv_obj <- CreateInfercnvObject(
                raw_counts_matrix = expr_matrix,
                annotations_file = annotations_file,
                delim = "\t",
                gene_order_file = gene_order_file,
                ref_group_names = ref_group_names,
                min_max_counts_per_cell = c(100, +Inf)  # Explicitly set cell filtering
            )
        })
    }, error = function(e) {
        self$say(paste("Error in CreateInfercnvObject:", e$message))
        stop(e)
    })
    
    # Extract and report final dimensions after InferCNV filtering
    final_cells <- ncol(infercnv_obj@expr.data)
    final_genes <- nrow(infercnv_obj@expr.data)
    
    cells_kept_pct <- round((final_cells / initial_cells) * 100, 1)
    genes_kept_pct <- round((final_genes / initial_genes) * 100, 1)
    
    self$say(paste("After filtering:", final_cells, "cells (", cells_kept_pct, "%) and", final_genes, "genes (", genes_kept_pct, "%) kept"))
    
    # Run infercnv analysis
    self$say("Running infercnv analysis...")
    
    # Set scientific notation options to prevent hclust errors
    old_scipen <- getOption("scipen")
    options(scipen = 100)
    on.exit(options(scipen = old_scipen), add = TRUE)
    
    suppressMessages({
        infercnv_obj <- infercnv::run(
            infercnv_obj,
            cutoff = cutoff,
            out_dir = out_dir,
            cluster_by_groups = cluster_by_groups,
            denoise = denoise,
            HMM = HMM,
            num_threads = num_threads,
            analysis_mode = analysis_mode
        )
    })
    
    # Store results in analysis
    self$analysis$infercnv <- infercnv_obj
    
    self$say(paste("InferCNV analysis completed. Results saved to:", out_dir))
    self$say("Results stored in self$analysis$infercnv")
    
    invisible(self)
})

ogel$set('public', 'plot_infercnv', function(
    plot_type = "observations",
    out_dir = NULL,
    output_filename = "infercnv_plot",
    output_format = "png",
    color_safe_pal = FALSE,
    x.center = 1,
    x.range = "auto",
    hclust_method = "ward.D2",
    custom_color_pal = NULL,
    plot_chr_scale = FALSE,
    chr_lengths = NULL,
    k_obs_groups = 1,
    contig_cex = 1
) {
    
    if (is.null(self$analysis$infercnv)) {
        stop("No InferCNV results found. Run infercnv() first.")
    }
    
    infercnv_obj <- self$analysis$infercnv
    
    if (is.null(out_dir)) {
        # Use default path since infercnv object structure varies
        out_dir <- file.path(self$path, "analysis", "infercnv")
    }
    
    # Ensure output directory exists
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    
    self$say("Plotting InferCNV results")
    
    # Plot function
    infercnv::plot_cnv(
        infercnv_obj,
        out_dir = out_dir,
        output_filename = output_filename,
        output_format = output_format,
        color_safe_pal = color_safe_pal,
        x.center = x.center,
        x.range = x.range,
        hclust_method = hclust_method,
        custom_color_pal = custom_color_pal,
        plot_chr_scale = plot_chr_scale,
        chr_lengths = chr_lengths,
        k_obs_groups = k_obs_groups,
        contig_cex = contig_cex
    )
    
    self$say(paste("Plot saved to:", file.path(out_dir, paste0(output_filename, ".", output_format))))
    
    invisible(self)
})

ogel$set('public', 'extract_cnv_regions', function(
    region_method = "HMM_CNV_predictions",
    min_cells_per_gene = 3,
    average_over_cells = TRUE,
    output_dir = NULL
) {
    
    if (is.null(self$analysis$infercnv)) {
        stop("No InferCNV results found. Run infercnv() first.")
    }
    
    infercnv_obj <- self$analysis$infercnv
    
    if (is.null(output_dir)) {
        # Use default path since infercnv object structure varies
        output_dir <- file.path(self$path, "analysis", "infercnv")
    }
    
    # Ensure output directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    self$say("Extracting CNV regions")
    
    # Extract CNV regions with error handling for different object structures
    tryCatch({
        if (region_method == "HMM_CNV_predictions") {
            if (!is.null(infercnv_obj@tumor_subclusters) && 
                "HMM_CNV_predictions" %in% names(infercnv_obj@tumor_subclusters)) {
                cnv_regions <- infercnv_obj@tumor_subclusters$HMM_CNV_predictions
            } else {
                self$say("Warning: HMM predictions not found, using expression matrix")
                cnv_regions <- infercnv_obj@expr.data
            }
        } else {
            cnv_regions <- infercnv_obj@expr.data
        }
    }, error = function(e) {
        self$say("Warning: Could not access expected slots, trying alternative access methods")
        # Fallback to main expression data
        if (!is.null(infercnv_obj@expr.data)) {
            cnv_regions <- infercnv_obj@expr.data
        } else {
            stop("Could not access CNV data from infercnv object")
        }
    })
    
    # Convert to data frame for easier handling
    cnv_df <- as.data.frame(as.matrix(cnv_regions))
    
    # Save CNV regions
    cnv_file <- file.path(output_dir, "cnv_regions.txt")
    write.table(cnv_df, file = cnv_file, sep = "\t", quote = FALSE, 
               row.names = TRUE, col.names = TRUE)
    
    # Note: Since we simplified storage, we can't store cnv_regions separately
    # Return the data frame for the user to handle
    self$say(paste("CNV regions extracted and saved to:", cnv_file))
    self$say(paste("CNV matrix dimensions:", nrow(cnv_df), "genes x", ncol(cnv_df), "cells"))
    
    return(cnv_df)
})

ogel$set('public', 'summarize_cnv', function(
    group_by = NULL,
    method = "mean"
) {
    
    if (is.null(self$analysis$infercnv)) {
        stop("No InferCNV results found. Run infercnv() first.")
    }
    
    # Always extract fresh CNV regions since we simplified storage
    self$say("Extracting CNV regions for summary...")
    cnv_df <- self$extract_cnv_regions()
    
    self$say("Summarizing CNV results")
    
    # Get metadata
    obj <- self$get_data()
    meta_data <- obj@meta.data
    
    # Match cells between CNV results and metadata
    common_cells <- intersect(colnames(cnv_df), rownames(meta_data))
    cnv_df <- cnv_df[, common_cells]
    meta_subset <- meta_data[common_cells, ]
    
    if (!is.null(group_by)) {
        if (!group_by %in% colnames(meta_subset)) {
            stop("group_by column '", group_by, "' not found in metadata")
        }
        
        groups <- meta_subset[[group_by]]
        unique_groups <- unique(groups)
        
        # Summarize by groups
        cnv_summary <- sapply(unique_groups, function(g) {
            group_cells <- common_cells[groups == g]
            group_cnv <- cnv_df[, group_cells, drop = FALSE]
            
            if (method == "mean") {
                rowMeans(group_cnv, na.rm = TRUE)
            } else if (method == "median") {
                apply(group_cnv, 1, median, na.rm = TRUE)
            } else {
                stop("method must be 'mean' or 'median'")
            }
        })
        
        cnv_summary <- as.data.frame(cnv_summary)
        colnames(cnv_summary) <- unique_groups
        
    } else {
        # Overall summary
        if (method == "mean") {
            cnv_summary <- rowMeans(cnv_df, na.rm = TRUE)
        } else if (method == "median") {
            cnv_summary <- apply(cnv_df, 1, median, na.rm = TRUE)
        }
        cnv_summary <- data.frame(cnv_score = cnv_summary)
    }
    
    # Note: Since we simplified storage, return summary instead of storing
    self$say("CNV summary completed")
    self$say(paste("Summary dimensions:", nrow(cnv_summary), "genes x", ncol(cnv_summary), "groups/conditions"))
    
    return(cnv_summary)
})

ogel$set('public', 'plot_cnv_downsample', function(
    n_cells = 3000,
    k_obs_groups = 3,
    output_filename = "infercnv_downsample",
    output_format = "pdf",
    out_dir = NULL,
    cluster_by_groups = TRUE,
    cluster_references = TRUE,
    plot_chr_scale = FALSE,
    x.center = 1,
    x.range = "auto",
    png_res = 300,
    useRaster = TRUE
) {
    
    if (is.null(self$analysis$infercnv)) {
        stop("No InferCNV results found. Run infercnv() first.")
    }
    
    if (is.null(out_dir)) {
        out_dir <- file.path(self$path, "analysis", "infercnv")
    }
    
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    
    self$say("Creating downsampled InferCNV plot")
    self$say(paste("Downsampling to", n_cells, "cells for visualization"))
    
    # Downsample the infercnv object
    infercnv_sample <- infercnv::sample_object(self$analysis$infercnv, n_cells = n_cells)
    
    # Plot the downsampled data
    infercnv::plot_cnv(
        infercnv_sample,
        k_obs_groups = k_obs_groups,
        cluster_by_groups = cluster_by_groups,
        cluster_references = cluster_references,
        plot_chr_scale = plot_chr_scale,
        out_dir = out_dir,
        x.center = x.center,
        x.range = x.range,
        output_filename = output_filename,
        output_format = output_format,
        png_res = png_res,
        useRaster = useRaster
    )
    
    self$say(paste("Downsampled plot saved as:", file.path(out_dir, paste0(output_filename, ".", output_format))))
    
    return(infercnv_sample)
})

ogel$set('public', 'infer_malignant_cells', function(
    top_percent = 0.05,
    score_threshold = 0.003,
    cor_threshold = 0.15,
    use_downsample = FALSE,
    n_sample_cells = 3000,
    plot_results = TRUE,
    plot_heatmap = TRUE,
    n_test_cells = 500,
    n_test_genes = 1000,
    out_dir = NULL
) {
    
    if (is.null(self$analysis$infercnv)) {
        stop("No InferCNV results found. Run infercnv() first.")
    }
    
    if (is.null(out_dir)) {
        out_dir <- file.path(self$path, "analysis", "infercnv")
    }
    
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    
    self$say("Inferring malignant cells using CNV patterns")
    self$say(paste("Parameters: top_percent =", top_percent, "score_threshold =", score_threshold, "cor_threshold =", cor_threshold))
    
    # Choose which infercnv object to use
    if (use_downsample) {
        self$say(paste("Using downsampled data with", n_sample_cells, "cells"))
        infercnv_obj <- infercnv::sample_object(self$analysis$infercnv, n_cells = n_sample_cells)
    } else {
        self$say("Using full dataset")
        infercnv_obj <- self$analysis$infercnv
    }
    
    # Get cell indices for tumor and normal cells
    if (length(infercnv_obj@observation_grouped_cell_indices) == 0) {
        stop("No observation groups found in InferCNV object")
    }
    if (length(infercnv_obj@reference_grouped_cell_indices) == 0) {
        stop("No reference groups found in InferCNV object")
    }
    
    # Get the first observation and reference groups
    tumor_cell_indices <- infercnv_obj@observation_grouped_cell_indices[[1]]
    normal_cell_indices <- infercnv_obj@reference_grouped_cell_indices[[1]]
    
    # Extract expression data
    obs <- infercnv_obj@expr.data[, tumor_cell_indices]
    ref <- infercnv_obj@expr.data[, normal_cell_indices]
    
    self$say(paste("Analyzing", ncol(obs), "tumor cells and", ncol(ref), "reference cells"))
    
    # Calculate CNV scores (deviation from 1.0)
    cnv_obs <- colMeans((obs - 1)^2)
    cnv_ref <- colMeans((ref - 1)^2)
    
    # Identify top tumor cells based on CNV score
    n_top_cells <- max(1, floor(length(cnv_obs) * top_percent))
    cell_top <- names(sort(cnv_obs, decreasing = TRUE))[1:n_top_cells]
    
    self$say(paste("Using top", n_top_cells, "cells (", top_percent*100, "%) as reference for correlation"))
    
    # Calculate reference profile from top cells
    if (length(cell_top) > 1) {
        cnv_top <- rowMeans(obs[, cell_top])
    } else {
        cnv_top <- obs[, cell_top]
    }
    
    # Calculate correlation with top tumor cells
    cor_obs <- apply(obs, 2, function(x) cor(x, cnv_top))
    cor_ref <- apply(ref, 2, function(x) cor(x, cnv_top))
    
    # Create results data frame
    cnv_results <- data.frame(
        score = c(cnv_obs, cnv_ref),
        cor = c(cor_obs, cor_ref),
        tissue = c(rep("tumor", length(cor_obs)), rep("normal", length(cor_ref))),
        stringsAsFactors = FALSE
    )
    rownames(cnv_results) <- c(names(cnv_obs), names(cnv_ref))
    
    # Classify malignant cells
    cnv_results$type <- "Not Malignant"
    cnv_results$type[cnv_results$cor > cor_threshold | cnv_results$score > score_threshold] <- "Malignant"
    cnv_results$type[cnv_results$tissue == "normal"] <- "Not Malignant"  # Force normal cells to be non-malignant
    
    # Summary statistics
    malignant_count <- sum(cnv_results$type == "Malignant" & cnv_results$tissue == "tumor")
    total_tumor <- sum(cnv_results$tissue == "tumor")
    malignant_pct <- round((malignant_count / total_tumor) * 100, 1)
    
    self$say(paste("Results:", malignant_count, "out of", total_tumor, "tumor cells classified as malignant (", malignant_pct, "%)"))
    
    # Plot results if requested
    if (plot_results) {
        suppressMessages({
            library(ggplot2)
            library(dplyr)
        })
        
        # Score vs correlation plot
        p1 <- cnv_results %>%
            ggplot(aes(score, cor)) +
            geom_point(aes(color = tissue), size = 0.5, alpha = 0.5) +
            geom_vline(xintercept = score_threshold, lty = "dashed") +
            geom_hline(yintercept = cor_threshold, lty = "dashed") +
            labs(title = "CNV Score vs Correlation with Top Tumor Cells",
                 x = "CNV Score", y = "Correlation") +
            theme_minimal()
        
        # Malignant classification plot
        p2 <- cnv_results %>%
            ggplot(aes(score, cor)) +
            geom_point(aes(color = type), size = 0.5, alpha = 0.5) +
            geom_vline(xintercept = score_threshold, lty = "dashed") +
            geom_hline(yintercept = cor_threshold, lty = "dashed") +
            scale_color_manual(values = c("Not Malignant" = "#8DD3C7", "Malignant" = "red")) +
            labs(title = "Malignant Cell Classification",
                 x = "CNV Score", y = "Correlation") +
            theme_minimal()
        
        # Save plots
        ggsave(file.path(out_dir, "malignant_cells_score_correlation.pdf"), p1, width = 6, height = 5)
        ggsave(file.path(out_dir, "malignant_cells_classification.pdf"), p2, width = 6, height = 5)
        
        self$say("Plots saved to output directory")
    }
    
    # Create heatmap if requested
    if (plot_heatmap && !is.null(infercnv_obj@gene_order)) {
        suppressMessages({
            library(ComplexHeatmap)
        })
        
        # Sample cells and genes for heatmap
        tumor_cells <- rownames(cnv_results)[cnv_results$tissue == "tumor"]
        set.seed(123)
        test_cells <- sample(tumor_cells, min(n_test_cells, length(tumor_cells)))
        set.seed(123)
        test_genes <- sample(rownames(infercnv_obj@gene_order), min(n_test_genes, nrow(infercnv_obj@gene_order)))
        
        # Create heatmap
        ht <- t(infercnv_obj@expr.data[test_genes, test_cells]) %>%
            Heatmap(
                name = "CNV",
                column_split = infercnv_obj@gene_order[test_genes, "chr"],
                cluster_columns = TRUE,
                cluster_column_slices = FALSE,
                show_column_names = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                cluster_row_slices = FALSE,
                right_annotation = rowAnnotation(
                    tissue = cnv_results[test_cells, "tissue"],
                    malignant = cnv_results[test_cells, "type"],
                    col = list(
                        tissue = c(tumor = "royalblue"),
                        malignant = structure(c("#8DD3C7", "red"),
                                            names = c("Not Malignant", "Malignant"))
                    )
                ),
                row_split = cnv_results[test_cells, "type"]
            )
        
        # Save heatmap
        pdf(file.path(out_dir, "malignant_cells_heatmap.pdf"), width = 12, height = 8)
        draw(ht)
        dev.off()
        
        self$say("Heatmap saved to output directory")
    }
    
    # Save results
    results_file <- file.path(out_dir, "malignant_cells_annotation.txt")
    write.table(cnv_results, results_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    self$say(paste("Results saved to:", results_file))
    
    return(cnv_results)
})

# 
# mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 
# genes <- getBM(
#   attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
#   mart = mouse_mart
# )
# 
# main_chrs <- c(1:19, "X", "Y")
# genes <- genes[genes$chromosome_name %in% main_chrs, ]
# genes$chromosome_name <- paste0("chr", genes$chromosome_name)
# gene_order <- genes[, c("external_gene_name", "chromosome_name", "start_position", "end_position")]
# 
# gene_order <- gene_order[!duplicated(gene_order$external_gene_name), ]
# gene_order <- gene_order[order(gene_order$chromosome_name, gene_order$start_position), ]
# write.table(gene_order, "mm10_gene_order.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
