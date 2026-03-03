# set plot functions here
ogel$set('public', 'plot_dim', function(group_by = self$celltype,label=T,label_insitu=T, reduction = "umap",palette = "seurat", theme = "theme_blank", ...,data_use = 'data') {
  #DimPlot(self[[data_use]], ...)
  CellDimPlot(self[[data_use]], group_by = group_by, reduction = reduction,theme = theme,label=label,label_insitu=label_insitu,...)
})

ogel$set('public','plot_features', function(features,order=T,label=T,...,data_use = 'data') {
  FeaturePlot(self[[data_use]],features = features,order = order,label=label,...)
})

ogel$set('public','plot_vln',function(...,data_use = 'data'){
  Seurat::VlnPlot(self[[data_use]],pt.size = 0,...)
})

ogel$set('public','plot_dot',function(features,reverse=F,size_range = c(0,6),group_by = self$celltype,features_df=NULL,group_df=NULL,data_use = 'data') {
  p_data <- self$get_features_expr_df(features=features,group_by = group_by,features_df=features_df,group_df=group_df,data_use = data_use)
  if(reverse){
    p <- 
      p_data%>%
      ggplot(aes(group_by,features))
  }else{
    p <- 
      p_data%>%
      ggplot(aes(features,group_by))
  }
  p +
    geom_point(aes(size = pct_exp,fill = avg_expr_scaled),
               color  = 'black',shape = 21) + 
    scale_size(range = size_range) + 
    scale_fill_gradientn(colors = c( '#2166ac','#67a9cf','#d1e5f0', "white", '#fddbc7', '#ef8a62', '#b2182b')) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill =NA,color = 'black'),
          axis.title = element_blank(),
          text = element_text(family = self$plot_config$fontfamily, size = self$plot_config$fontsize)) +
    labs(fill = 'Scaled expression',size = "Percent expressed")
})

ogel$set('public','plot_feature_density',function(features,data_use = 'data',...) {
  Nebulosa::plot_density(self[[data_use]], features,...)
})

# Input: diff_table - data frame with differential analysis results
# Output: ggplot2 volcano plot object
ogel$set('public', 'plot_jj_volcano', function(diff_table,
                                                feature_col = 'feature',
                                                cluster_col = 'group',
                                                logFC_col = 'logFC',
                                                pvalue_col = 'pval',
                                                padj_col = 'padj',
                                                base_size = 10,
                                                topGeneN = 3,
                                                plot_name = 'jj_volcano',
                                                log2FC.cutoff = 0.25,
                                                pvalue.cutoff = 0.05,
                                                adjustP.cutoff = 0.01,
                                                col.type = 'updown',
                                                back.col = 'grey93',
                                                pSize = 0.75,
                                                aesCol = c('#0099CC', '#CC3333'),
                                                legend.position = 'top',
                                                flip = FALSE,
                                                celltypeSize = 3,
                                                cluster_cols = NULL,
                                                group_labels = c('Control up', 'LEB up'),
                                                ylab_text = 'Average log2FoldChange',
                                                label_size = 3,
                                                ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) stop('ggplot2 required')
  if (!requireNamespace('ggrepel', quietly = TRUE)) stop('ggrepel required')
  if (is.null(diff_table) || nrow(diff_table) == 0) return(invisible(NULL))

  df <- data.frame(
    gene = diff_table[[feature_col]],
    cluster = diff_table[[cluster_col]],
    avg_log2FC = diff_table[[logFC_col]],
    p_val = diff_table[[pvalue_col]],
    p_val_adj = diff_table[[padj_col]],
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0) return(invisible(NULL))

  # Filter significant rows else relax
  sig <- df[abs(df$avg_log2FC) >= log2FC.cutoff & df$p_val < pvalue.cutoff, , drop = FALSE]
  if (nrow(sig) == 0) {
    sig <- df
    log2FC.cutoff <- min(0.1, stats::quantile(abs(df$avg_log2FC), 0.75, na.rm = TRUE))
    pvalue.cutoff <- max(0.1, stats::quantile(df$p_val, 0.75, na.rm = TRUE))
  }
  sig$type <- ifelse(sig$avg_log2FC >= log2FC.cutoff, group_labels[2], group_labels[1])
  sig$type2 <- ifelse(sig$p_val_adj < adjustP.cutoff,
                      paste('adjust Pvalue <', adjustP.cutoff),
                      paste('adjust Pvalue >=', adjustP.cutoff))

  # Background extent per cluster
  available_clusters <- unique(sig$cluster)
  back.data <- do.call(rbind, lapply(available_clusters, function(x) {
    tmp <- sig[sig$cluster == x, , drop = FALSE]
    if (nrow(tmp) > 0) {
      data.frame(cluster = x,
                 min = min(tmp$avg_log2FC, na.rm = TRUE) - 0.2,
                 max = max(tmp$avg_log2FC, na.rm = TRUE) + 0.2)
    } else {
      data.frame(cluster = x, min = -0.5, max = 0.5)
    }
  }))

  # Select top markers per cluster
  safe_slice <- function(d, n, decreasing = TRUE) {
    if (nrow(d) == 0) return(d)
    ord <- if (decreasing) order(-d$avg_log2FC) else order(d$avg_log2FC)
    head(d[ord, , drop = FALSE], n)
  }
  top.marker <- do.call(rbind, lapply(split(sig, sig$cluster), function(d) {
    rbind(safe_slice(d, topGeneN, TRUE), safe_slice(d, topGeneN, FALSE))
  }))

  # Colors for clusters (tiles)
  cluster_levels <- unique(as.character(sig$cluster))
  if (is.null(cluster_cols)) {
    cluster_cols <- grDevices::rainbow(length(cluster_levels))
    names(cluster_cols) <- cluster_levels
  } else {
    if (!all(cluster_levels %in% names(cluster_cols))) {
      missing <- setdiff(cluster_levels, names(cluster_cols))
      warning(paste("Missing colors for:", paste(missing, collapse = ", ")))
      cluster_cols <- c(cluster_cols, setNames(grDevices::rainbow(length(missing)), missing))
    }
    cluster_cols <- cluster_cols[cluster_levels]
  }

  # Build plot
  p1 <- ggplot2::ggplot(sig, ggplot2::aes(x = cluster, y = avg_log2FC)) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = min), fill = back.col) +
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, y = max), fill = back.col)

  if (col.type == 'updown') {
    color_values <- setNames(aesCol, group_labels)
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize) +
      ggplot2::scale_color_manual(values = color_values)
  } else if (col.type == 'adjustP') {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2), size = pSize) +
      ggplot2::scale_color_manual(values = c(aesCol[2], aesCol[1]))
  } else {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize)
  }

  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank()) +
    ggplot2::xlab('Clusters') +
    ggplot2::ylab(ylab_text) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, fill = cluster), color = 'black',
                       height = log2FC.cutoff * 2, alpha = 0.3, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cluster_cols)

  if (!is.null(top.marker) && nrow(top.marker) > 0) {
    p4 <- p4 + ggrepel::geom_text_repel(data = top.marker,
                                        ggplot2::aes(x = cluster, y = avg_log2FC, label = gene),
                                        size = label_size,
                                        max.overlaps = 50, ...)
  }

  p5 <- if (isTRUE(flip)) {
    p4 + ggplot2::scale_y_continuous(n.breaks = 6) +
      ggplot2::geom_label(ggplot2::aes(x = cluster, y = 0, label = cluster), size = celltypeSize) +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::coord_flip()
  } else {
    p4 + ggplot2::geom_label(ggplot2::aes(x = cluster, y = 0, label = cluster), size = celltypeSize) +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  }

  if (self$auto_save) self$save_fig(p5, plot_name, width = 8, height = 6)
  p5
})