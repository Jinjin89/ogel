ogel$set('public','diff_wilcoxauc',function(group_by,..., assay_use='RNA',force=F,results_name = NULL,return_df=F,data_use = 'data'){
  if(!'wilcoxauc' %in% names(self$analysis)){
    self$analysis$wilcoxauc <- list()
  }
  
  if(is.null(results_name)){
    results_name <- group_by
  }
  raw_levels <- levels(factor(self$get_data(data_use)@meta.data[[group_by]]))
  if((!results_name %in% names(self$analysis$wilcoxauc)) || force){
    message('runing wilcoxauc for: ',group_by)
    df <- presto::wilcoxauc(
      self$get_data(data_use),
      seurat_assay = assay_use,
      ...,
      group_by = group_by)  %>% 
      dplyr::mutate(pct_change = pct_in - pct_out) %>% 
      dplyr::mutate(group = factor(group, levels = raw_levels))
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
  warning('TODO: suppport multiple group heatmap')
  
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
.enrich_gsea <- \(enrich_data,input_db_list,pvalueCutoff = 1,eps = 0){
  db_names <- names(input_db_list)
  sifn(len(db_names) > 0, 'input db should be with names')
  enrich_res <- list()
  for(each_db_name in db_names){
    cat("GSEA: the db is: ",each_db_name,'\n')
    enrich_res[[each_db_name]] <- clusterProfiler::GSEA(geneList = enrich_data, TERM2GENE = input_db_list[[each_db_name]], pvalueCutoff = pvalueCutoff, eps = 0)  
  }
  return(enrich_res)
}

.enrich_ora <- \(enrich_data,input_db_list,pvalueCutoff = 1){
  db_names <- names(input_db_list)
  sifn(len(db_names) > 0, 'input db should be with names')
  enrich_res <- list()
  for(each_db_name in db_names){
    cat("ORA: the db is: ",each_db_name,'\n')
    enrich_res[[each_db_name]] <- clusterProfiler::enricher(
      gene = enrich_data,
      TERM2GENE = input_db_list[[each_db_name]],
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = pvalueCutoff )
    
  }
  return(enrich_res)
}

ogel$set('public','diff_enrich_gsea',function(
    results_name,enrich_res_name = results_name,
    group_show= NULL,db_name = c('hallmark','kegg','go'),logFC = 'logFC',
    feature = 'feature',pvalueCutoff = 1,return_enrich = F,force = F){
  stopifnot(all(db_name %in% names(self$db)))
  db <- self$db[db_name]
  if(!'enrich' %in% names(self$analysis)){
    self$analysis$enrich <- list()
  }
  message('enrichment analysis for: ', results_name, ' using db: ', db_name)
  if(is.null(self$analysis$wilcoxauc[[results_name]])){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }
  
  diff_res <- self$analysis$wilcoxauc[[results_name]]
  all_groups_found <- diff_res$group %>% unique() %>% as.character()
  enrich_res <- list()
  if((!enrich_res_name %in% names(self$analysis$enrich)) || force){
    if(len(group_show) == 0){
      group_show <- all_groups_found
    }else{
      group_show <- intersect(group_show,all_groups_found)
    }
    sifn(len(group_show) > 0,'should have at least group to enrich')
    # for loop for all
    for(each_group in group_show){
      self$sep_line()
      cat("Enrich for: ",each_group,'\n')
      # 1) get the genes for supplied groups
      input_degs <- diff_res %>% dplyr::filter(group == each_group)
      # 3) genes logFC
      input_genes = dplyr::arrange(input_degs, dplyr::desc(!!as.name(logFC)))
      genes = input_genes[[logFC]]
      names(genes) = input_genes[[feature]]
      enrich_res[[each_group]] = .enrich_gsea(enrich_data = genes,input_db_list = db,pvalueCutoff = pvalueCutoff,eps = 0)
    }
    self$analysis$enrich[[enrich_res_name]] <- enrich_res
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

ogel$set('public','diff_enrich_ora',function(
    results_name,enrich_res_name = results_name, top = 100,diff_padj_cutoff=0.01,
    group_show= NULL,db_name = c('hallmark','kegg','go'),logFC = 'logFC',
    feature = 'feature',pvalueCutoff = 1,return_enrich = F,force = F){
  stopifnot(all(db_name %in% names(self$db)))
  db <- self$db[db_name]
  if(!'enrich' %in% names(self$analysis)){
    self$analysis$enrich <- list()
  }
  message('enrichment analysis for: ', results_name, ' using db: ', db_name)
  if(is.null(self$analysis$wilcoxauc[[results_name]])){
    stop('wilcoxauc results not found for group: ', results_name, '. Run diff_wilcoxauc first.')
  }
  
  diff_res <- self$analysis$wilcoxauc[[results_name]]
  all_groups_found <- diff_res$group %>% unique() %>% as.character()
  enrich_res <- list()
  if((!enrich_res_name %in% names(self$analysis$enrich)) || force){
    if(len(group_show) == 0){
      group_show <- all_groups_found
    }else{
      group_show <- intersect(group_show,all_groups_found)
    }
    sifn(len(group_show) > 0,'should have at least group to enrich')
    
    # for loop for all
    for(each_group in group_show){
      self$sep_line()
      cat("Enrich for: ",each_group,'\n')
      # 1) get the genes for supplied groups
      input_degs <- diff_res %>% dplyr::filter(group == each_group) %>% 
        dplyr::filter(padj < diff_padj_cutoff)
      if(top >0){
        input_degs <- 
          input_degs %>% 
          dplyr::arrange(dplyr::desc(!!as.name(logFC))) %>% 
          head(top) %>% 
          dplyr::pull(feature)
      }else{
        input_degs <- 
          input_degs %>% 
          dplyr::arrange(!!as.name(logFC)) %>% 
          head(abs(top)) %>% 
          dplyr::pull(feature)
      }
      # 3) genes logFC
      enrich_res[[each_group]] = .enrich_ora(enrich_data = input_degs,input_db_list = db,pvalueCutoff = pvalueCutoff)
    }
    self$analysis$enrich[[enrich_res_name]] <- enrich_res
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

ogel$set('public','diff_top_features',function(results_name,top = 5,
  by_what = 'pct_change',pvalue_use = 'padj',
  min_pct_in = 20,
  filter_weird = F,
  pvalue_cutoff = 0.01){
  diff_table_by_wilcoxauc <- self$analysis$wilcoxauc[[results_name]] %>% 
    dplyr::filter(pct_in > min_pct_in)
  # This function return a list of genes
  if(filter_weird){
    cat("WARNING: `filter_weird` is set to TRUE, This may remove some genes that you are interested!\n")
    diff_table_by_wilcoxauc <- 
      diff_table_by_wilcoxauc %>% 
      dplyr::filter(!stringr::str_detect(feature,'^\\d')) %>% 
      dplyr::filter(!stringr::str_detect(feature,'\\.\\d$'))
  }
  diff_table_by_wilcoxauc <- 
    diff_table_by_wilcoxauc %>% 
    dplyr::filter(.[[pvalue_use]] < pvalue_cutoff) %>% 
    dplyr::group_by(group) %>% 
    dplyr::arrange(dplyr::desc(!!as.name(by_what))) %>% 
    do(head(.,abs(top))) %>% 
    ungroup() %>% 
    dplyr::filter(!duplicated(feature))
  diff_table_by_wilcoxauc
}

)
ogel$set("public",'diff_plot_hm', function(features, group_by = self$celltype,data_use = 'data') {
  ComplexHeatmap::ht_opt(message = F)
  # downsample_for ploting
  dta <- self$get_data(data_use)
  Idents(dta) <- group_by
  set.seed(1)

  dta <- subset(x = dta, downsample = 500)
  
  unique_features <- unique(features) %>% dplyr::intersect(rownames(dta))
  anno_top <- dta@meta.data[,group_by,drop=F]
  all_val <- sort(unique(as.character(anno_top[[1]])))
  
  col_list <- list(g = structure(.scPalette1(length(all_val)),names = all_val))
  names(col_list) <- group_by
  
  # replot using complexhemap
  temp_expr <- GetAssayData(dta,layer='data')[unique_features,]
  temp_expr <- t(apply(as.matrix(temp_expr),1,scale))
  #temp_expr <- as.matrix(GetAssayData(dta,layer='data')[unique_features,])
  temp_expr <- na.omit(temp_expr)
  colnames(temp_expr) <- colnames(dta)
  
  # Determine which rows to label with anno_mark
  if (nrow(temp_expr) > 35) {
    set.seed(12)
    mark_indices <- sample(seq_len(nrow(temp_expr)), 35)
  } else {
    mark_indices <- seq_len(nrow(temp_expr))
  }
  suppressWarnings(
    {
      hm <- 
        ComplexHeatmap::Heatmap(
          temp_expr,
          # row
          row_names_side = 'left',
          show_row_names = FALSE, # Hide default row names
          column_split = anno_top[[group_by]],
          show_row_dend = F,
          cluster_rows = F,
          
          # Use row annotation mark for showing selected row labels
          left_annotation = ComplexHeatmap::rowAnnotation(
            mark = ComplexHeatmap::anno_mark(
              at = mark_indices,
              labels = rownames(temp_expr)[mark_indices], 
              side = "left", 
              labels_gp = grid::gpar(fontsize = 10),
              link_width = grid::unit(5, "mm")
            )
          ),
          
          # col
          show_column_names = F,
          show_column_dend = F,
          column_title = ' ',
          cluster_column_slices = F,
          top_annotation = HeatmapAnnotation(
            df = anno_top,
            col = col_list,
            show_annotation_name = F),
          
          # main params
          name = 'Z-score',
          border = 'black',
          col = circlize::colorRamp2(
            breaks = seq(-3,3,length.out = 5),
            colors = c( '#2166ac','#67a9cf', "white", '#ef8a62', '#b2182b'))
        )
    }
  )
  cat("recommended save plot size is(width X height): 8 X 6\n")
  hm
}
)


ogel$set('public','diff_plot_dot',function(
  features,
  group_by = self$celltype,
  left_annotation_color = NULL,
  dot_size_scale = 5,
  scale_min = -1,
  scale_max = 1,
  scale_colors = c("#2166ac", "white", "#b2182b"),
  # legend
  fill_label = 'Average Expression',
  size_label = 'Percent Expressed',
  data_use = 'data'
) {
  # Load required packages
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  
  # 1) Extract data using DotPlot (same as original function)
  obj <- self$get_data(data_use)
  data_get <- Seurat::DotPlot(obj, features = features, group.by = group_by)
  df <- data_get$data
  
  # Make sure gene ordering is preserved
  df$features.plot <- factor(df$features.plot, levels = rev(unique(df$features.plot)))
  
  # Cap expression values for better visualization
  df$avg.exp.scaled[df$avg.exp.scaled > scale_max] <- scale_max
  df$avg.exp.scaled[df$avg.exp.scaled < scale_min] <- scale_min
  
  # Generate colors for clusters
  if (is.null(left_annotation_color)) {
    n_clusters <- length(unique(df$id))
    cluster_colors <- scales::hue_pal()(n_clusters)
    names(cluster_colors) <- unique(df$id)
  } else {
    cluster_colors <- left_annotation_color
    names(cluster_colors) <- unique(df$id)
  }
  
  # 2) Create the ggplot2 dotplot
  p <- ggplot(df, aes(y = id, x = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_size_continuous(
      range = c(0, dot_size_scale),
      name = "Percent Expressed"
    ) +
    scale_color_gradient2(
      low = scale_colors[1],
      mid = scale_colors[2],
      high = scale_colors[3],
      midpoint = 0,
      name = "Average Expression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      panel.border = element_rect(fill = NA, color = 'black'),
      axis.title.y = element_blank(),
      legend.position = "top"
    ) +
    labs(fill = fill_label,size = size_label)
  
  # which size to use
  # 
  plot_width <- max(8, 2.5 + 0.25 * length(unique(features)))
  plot_height <- max(4, 2 + 0.4 * length(unique(obj@meta.data[[group_by]])))
  cat("infer save size using:\n-> max(8, 2.5 + 0.25 * length(unique(features)))\n-> max(4, 2 + 0.4 * length(unique(obj@meta.data[[group_by]])))\n")
  cat("recommended save plot size is(width X height): ", plot_width, 'X' ,plot_height,'\n')
  
  return(p)
}
)


.scPalette1 <- function (n) {
  colorSpace <- c("#377EB8", "#4DAF4A", "#984EA3", 
                  "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                  "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                  "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                  "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

