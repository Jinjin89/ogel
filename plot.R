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
          axis.title = element_blank()) +
    labs(fill = 'Scaled expression',size = "Percent expressed")
})

ogel$set('public','plot_feature_density',function(features,data_use = 'data',...) {
  Nebulosa::plot_density(self[[data_use]], features,...)
})