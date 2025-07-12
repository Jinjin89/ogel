ogel$set('public','diff_wilcoxauc',function(group_by,..., assay_use='RNA',force=F,results_name = NULL,return_df=F,data_use = 'data'){
  if(!'wilcoxauc' %in% names(self$analysis)){
    self$analysis$wilcoxauc <- list()
  }
  if(is.null(results_name)){
    results_name <- group_by
  }

  if((!results_name %in% names(self$analysis$wilcoxauc)) || force){
    message('runing wilcoxauc for: ',group_by)
    df <- presto::wilcoxauc(
      self[[data_use]],
      seurat_assay = assay_use,
      ...,
      group_by = group_by)  %>% 
      dplyr::mutate(pct_change = pct_in - pct_out) %>% 
      dplyr::mutate(group = forcats::fct_inorder(group))
      self$analysis$wilcoxauc[[results_name]] <- df
    message('wilcoxauc results saved to: ',results_name)
  }else{
    message('wilcoxauc results found! not running')
  }
  if(return_df){
    df <- self$analysis$wilcoxauc[[results_name]]
    return(df)
  }else{
    invisible(self)
  }
  
})

ogel$set('public','diff_plot_volcano',function(results_name,group_show,logFC_cutoff = 0.25, pval_cutoff = 0.05, 
                                      top_n = 8, feature_include = NULL,
                                      direction = c("both", "y", "x"), label_size = 3,
                                      maxup_add = 1, maxlower_add = 1){
  message('volcano plot from wilcoxauc for: ', results_name)
  
  # Check if wilcoxauc results exist
  if(!'wilcoxauc' %in% names(self$analysis) || 
     !results_name %in% names(self$analysis$wilcoxauc)){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }

  input_degs <- self$analysis$wilcoxauc[[results_name]] %>%
    dplyr::filter(group == group_show)
  
  # Add p_sig column if not present
  if(!'p_sig' %in% colnames(input_degs)){
    input_degs <- input_degs %>%
      dplyr::mutate(
        p_sig = dplyr::case_when(
          padj < pval_cutoff & logFC > logFC_cutoff ~ "Up",
          padj < pval_cutoff & logFC < -logFC_cutoff ~ "Down",
          TRUE ~ "N.S."
        )
      )
  }
  
  stopifnot('p_sig' %in% colnames(input_degs))
  
  # Replace zero p-values
  padj_replace <- input_degs %>% dplyr::filter(padj > 0) %>% dplyr::pull(padj) %>% min
  input_degs$padj <- ifelse(input_degs$padj == 0, padj_replace, input_degs$padj)
  
  # Get top markers
  markers_top <- input_degs %>% 
    dplyr::filter(p_sig == "Up") %>% 
    dplyr::arrange(dplyr::desc(abs(logFC))) %>% 
    dplyr::mutate(i = 1:dplyr::n()) %>% 
    dplyr::mutate(labels = ifelse(i <= top_n | feature %in% feature_include, feature, NA))
    
  markers_down <- input_degs %>% 
    dplyr::filter(p_sig == "Down") %>% 
    dplyr::arrange(dplyr::desc(abs(logFC))) %>% 
    dplyr::mutate(i = 1:dplyr::n()) %>% 
    dplyr::mutate(labels = ifelse(i <= top_n | feature %in% feature_include, feature, NA))
  
  # Create volcano plot
  p <- input_degs %>% 
    ggplot2::ggplot(ggplot2::aes(x = logFC, y = -log10(padj))) + 
    ggplot2::geom_point(size = 2, alpha = 0.2, ggplot2::aes(color = p_sig)) + 
    ggplot2::geom_vline(xintercept = abs(logFC_cutoff), color = "gray", lty = "dashed") +
    ggplot2::geom_vline(xintercept = -abs(logFC_cutoff), color = "gray", lty = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(pval_cutoff), color = "gray", lty = "dashed") +
    ggplot2::scale_color_manual(values = c('Down' = 'blue', 'N.S.' = 'gray', 'Up' = 'red')) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +
    ggrepel::geom_label_repel(
      data = markers_top %>% dplyr::filter(logFC > 0) %>% tidyr::drop_na(labels),
      ggplot2::aes(label = labels),
      segment.color = 'gray',
      color = 'red',
      size = label_size,
      max.overlaps = Inf,
      direction = direction[1],
      vjust = 0
    ) +
    ggrepel::geom_label_repel(
      data = markers_down %>% tidyr::drop_na(labels),
      ggplot2::aes(label = labels),
      segment.color = 'gray',
      color = 'blue',
      size = label_size,
      max.overlaps = Inf,
      direction = direction[1],
      vjust = 0
    )
  
  return(p)
})

ogel$set('public','diff_plotheatmap',function(results_name,group_by,group_show,..., logFC_cutoff = 0.25, pval_cutoff = 0.05,
                                       downsample = NULL, top_n = 8, feature_include = NULL, 
                                       plot_genes = 50, data_use = 'data'){
  message('heatmap plot from wilcoxauc for: ', results_name)
  
  # Check if wilcoxauc results exist
  if(!'wilcoxauc' %in% names(self$analysis) || 
     !results_name %in% names(self$analysis$wilcoxauc)){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }
  
  input_degs <- self$analysis$wilcoxauc[[results_name]] %>%
    dplyr::filter(group == group_show)
  
  # Add p_sig column if not present
  if(!'p_sig' %in% colnames(input_degs)){
    input_degs <- input_degs %>%
      dplyr::mutate(
        p_sig = dplyr::case_when(
          padj < pval_cutoff & logFC > logFC_cutoff ~ "Up",
          padj < pval_cutoff & logFC < -logFC_cutoff ~ "Down",
          TRUE ~ "N.S."
        )
      )
  }
  
  stopifnot('p_sig' %in% colnames(input_degs))
  
  # Function to select markers
  fun_select_markers <- function(input_df, plot_genes){
    input_df$pct_change <- input_df$pct_in - input_df$pct_out
    input_df <- input_df %>% dplyr::filter(p_sig != "N.S.")
    
    genes_up <- input_df %>% 
      dplyr::arrange(dplyr::desc(pct_change)) %>% 
      utils::head(plot_genes) %>% 
      dplyr::filter(logFC > 0) %>% 
      dplyr::pull(feature)
    
    genes_down <- input_df %>% 
      dplyr::arrange(pct_change) %>% 
      dplyr::filter(logFC < 0) %>% 
      utils::head(plot_genes) %>% 
      dplyr::pull(feature)
    
    c(genes_up, genes_down) %>% unique()
  }
  
  # Get markers for heatmap
  markers_find <- fun_select_markers(input_degs, plot_genes = plot_genes)
  markers_for_heatmap <- c(feature_include, markers_find) %>% unique()
  
  # Get object for plotting
  obj_data <- self[[data_use]]
  obj_meta <- obj_data@meta.data
  
  # Downsample if specified
  if(!is.null(downsample) && length(downsample) > 0){
    message('*   downsampling: ', downsample)
    # Simple downsampling by group
    if(group_by %in% colnames(obj_meta)){
      sampled_cells <- obj_meta %>%
        dplyr::group_by(!!rlang::sym(group_by)) %>%
        dplyr::slice_sample(n = min(downsample, dplyr::n())) %>%
        dplyr::pull(rownames(.))
      obj_data <- obj_data[, sampled_cells]
      obj_meta <- obj_meta[sampled_cells, ]
    }
  }
  
  # Get expression data
  obj_expr <- GetAssayData(obj_data[intersect(markers_for_heatmap, rownames(obj_data)), ]) %>% 
    as.matrix() %>% 
    t() %>% 
    scale() %>% 
    t()
  
  # Remove genes with all NAs
  na_index <- apply(!is.na(obj_expr), 1, all)
  obj_expr <- obj_expr[unname(na_index), ]
  
  # Create annotation
  top_anno_df <- obj_meta %>% 
    dplyr::select(!!rlang::sym(group_by))
  
  # Create row annotation for top markers
  show_marker_index <- c(1:min(top_n, nrow(obj_expr)), 
                        max(1, nrow(obj_expr) - top_n + 1):nrow(obj_expr))
  show_marker_index <- unique(show_marker_index)
  
  right_annotation <- data.frame(
    feature = rownames(obj_expr),
    i = 1:nrow(obj_expr)
  ) %>% 
    dplyr::mutate(labels = ifelse(i %in% show_marker_index, feature, NA)) %>% 
    tidyr::drop_na(labels)
  
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    obj_expr,
    #col = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c('blue', 'white', 'red')),
    col = circlize::colorRamp2(breaks = c(-3,-2,-1, 0,1, 2,3), colors = c( '#2166ac','#67a9cf','#d1e5f0', "white", '#fddbc7', '#ef8a62', '#b2182b')),
    #col = c( '#2166ac','#67a9cf','#d1e5f0', "white", '#fddbc7', '#ef8a62', '#b2182b'),
    name = 'z-score',
    # Row settings
    show_row_dend = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    right_annotation = ComplexHeatmap::rowAnnotation(
      f = ComplexHeatmap::anno_mark(
        at = right_annotation$i,
        labels = right_annotation$feature
      )
    ),
    row_km = 2,
    # Column settings
    show_column_names = FALSE,
    show_column_dend = TRUE,
    cluster_column_slices = FALSE,
    cluster_columns = TRUE,
    column_split = top_anno_df[[group_by]],
    ...
  )
  
  return(ht)
})


# enrich step
ogel$set('public','diff_enrich_gsea',function(results_name,group_show,db_name,logFC = 'logFC',feature = 'feature',pvalueCutoff = 1,return_enrich = F,force = F){
  stopifnot(db_name %in% names(self$db))
  db <- self$db[[db_name]]
  if(!'enrich' %in% names(self$analysis)){
    self$analysis$enrich <- list()
  }
  message('enrichment analysis for: ', results_name, ' using db: ', db_name)
  if(is.null(self$analysis$wilcoxauc[[results_name]])){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }

  enrich_res_name <- paste0(results_name,'_',db_name)
  if((!enrich_res_name %in% names(self$analysis$enrich)) || force){
    # 1) get the genes for supplied groups
    input_degs <- self$analysis$wilcoxauc[[results_name]] %>%
      dplyr::filter(group == group_show)

    # 2) term 2 gene
      term_2_genes = dplyr::select(db, dplyr::all_of(c("term","gene")))

    # 3) genes logFC
      input_genes = dplyr::arrange(input_degs, dplyr::desc(!!as.name(logFC)))
      genes = input_genes[[logFC]]
      names(genes) = input_genes[[feature]]

      enrich_results = clusterProfiler::GSEA(geneList = genes, TERM2GENE = term_2_genes, pvalueCutoff = pvalueCutoff, eps = 0)
      self$analysis$enrich[[enrich_res_name]] <- enrich_results
      message('enrichment results saved to: ', enrich_res_name)
    }else{
      message('enrichment results found! not running')
    }
    if(return_enrich){
      return(self$analysis$enrich[[enrich_res_name]])
    }else{
      invisible(self)
    }
})

ogel$set('public','diff_enrich_gsea_patch',function(results_name,group_show,db_names = c('go','kegg','hallmark'),logFC = 'logFC',feature = 'feature',pvalueCutoff = 1,force = F){
  stopifnot(all(db_names %in% names(self$db)))
    if(!'enrich' %in% names(self$analysis)){
    self$analysis$enrich <- list()
  }
  if(is.null(self$analysis$wilcoxauc[[results_name]])){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }
  if(results_name %in% names(self$analysis$enrich) && !force){
    message('enrichment results found and force is FALSE! not running')
    return(self)
  }
  # 1) get the data for enriment
  input_degs <- self$analysis$wilcoxauc[[results_name]] %>%
    dplyr::filter(group == group_show) %>%
    dplyr::filter(!is.na(!!as.name(logFC))) %>%
    dplyr::filter(!is.na(!!as.name(feature))) %>%
    dplyr::filter(!!as.name(logFC) != 0)
    input_genes = dplyr::arrange(input_degs, dplyr::desc(!!as.name(logFC)))
    genes = input_genes[[logFC]]
    names(genes) = input_genes[[feature]]

  # 2) the enrich list with results_name
  self$analysis$enrich[[results_name]] <- list()
  self$analysis$enrich[[results_name]]$degs <- input_degs

  # 3) get the db for enriment
  for(db_name in db_names){
    message('enrichment analysis for: ', results_name, ' using db: ', db_name)
    term_2_genes = dplyr::select(self$db[[db_name]], dplyr::all_of(c("term","gene"))) %>%
      dplyr::filter(gene %in% names(genes))
    suppressWarnings({
    suppressMessages(
      {
        enrich_results = clusterProfiler::GSEA(geneList = genes, TERM2GENE = term_2_genes, pvalueCutoff = pvalueCutoff, eps = 0)
      })
    })
    self$analysis$enrich[[results_name]][[db_name]] <- enrich_results
  }# end for db_names loop
  invisible(self)
})
