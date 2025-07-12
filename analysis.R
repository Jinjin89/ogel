
ogel$set('public','stat_metadata_counts',function(
    by,celltype = self$celltype,
    by_remove = NULL, celltype_remove = NULL,
    by_features = NULL,celltype_features=NULL,
    fun_agg = 'mean',
    data_use = 'data'){
  stopifnot(len(by) == 1)
  
  meta.data <- self$get_data(data_use)@meta.data
  unique_by_total <- unique(meta.data[[by]])
  unique_celltype_total <- unique(meta.data[[celltype]])
  
  if(len(by_remove) >0 ){
    meta.data <- 
      meta.data %>% 
      dplyr::filter(!.[[by]] %in% by_remove)  
    
  }
  
  if(len(celltype_remove) >0 ){
    meta.data <- 
      meta.data %>% 
      dplyr::filter(!.[[celltype]] %in% celltype_remove)
  }

  cat("keep ",len(unique(meta.data[[by]])),'/',len(unique_by_total),' for by\n')
  cat("keep ",len( unique(meta.data[[celltype]])),'/',len(unique_celltype_total),' for celltypes \n')
  
  fun_collect_anno <- \(input_df,id,vars){
    stopifnot(is.data.frame(input_df))
    stopifnot(all(vars %in% colnames(input_df)))
    input_df <- input_df %>% dplyr::select(all_of(c(id,vars)))
    fun_agg_wrapper <- \(x){
      if(is.numeric(x)){
        return(do.call(fun_agg,list(x)))
      }else{
        x = unique(x)
        stopifnot(len(x) == 1)
        return(x)
      }
    }
    
    input_df %>% 
      dplyr::reframe(across(all_of(vars),fun_agg_wrapper),.by = all_of(id)) %>% 
      tibble::column_to_rownames(id)
  }
  
  if(len(celltype_features) >0){
    anno_celltype <- meta.data %>% 
      fun_collect_anno(celltype,celltype_features)
      
  }
  if(len(by_features) >0){
    anno_by <- meta.data %>% 
      fun_collect_anno(by,by_features)
  }
  
  # anno_by <- meta.data
  all_features_combine <- c(by,celltype,by_features,celltype_features)
  features_no_celltype <- all_features_combine %>% dplyr::setdiff(celltype)
  
  # 1) 
  stat_df <- 
    meta.data %>% 
    dplyr::count(across(all_of(c(by,celltype)))) %>% 
    tidyr::pivot_wider(id_cols = all_of(by),names_from = all_of(celltype),values_from = n,values_fill = 0) %>% 
    tidyr::pivot_longer(-all_of(by),names_to = 'celltype',values_to = 'n')
  
  for(i in seq_along(by_features)){
    stat_df[[by_features[i]]] <- anno_by[as.character(stat_df[[by]]),by_features[i]]
  }
  
  for(i in seq_along(celltype_features)){
    stat_df[[celltype_features[i]]] <- anno_celltype[as.character(stat_df[[celltype]]),celltype_features[i]]
  }
  
  stat_df <- 
    stat_df %>% 
    dplyr::mutate(pct = n/sum(n) * 100,.by = all_of(by))
  
  # 2)
  stat_mat <- stat_df %>% 
    tidyr::pivot_wider(names_from = all_of(celltype),id_cols = all_of(by),values_from = pct) %>% 
    tibble::column_to_rownames(by) %>% 
    as.matrix()
  
  res <- list(df = stat_df,mat = stat_mat)
  if(len(by_features) > 0){
    res[['anno_by']] <- anno_by[rownames(stat_mat),,drop=F]
  }
  if(len(celltype_features) > 0){
    res[['anno_celltype']] <- anno_celltype[colnames(stat_mat),,drop=F]
  }
  return(res)
})


ogel$set('public','plot_stat_metadata_counts',\(
  res,
  col = c( '#2166ac','#67a9cf','#d1e5f0','#fddbc7', '#ef8a62', '#b2182b'),
  right_annotation = NULL,
  row_bar_color = NULL,
  row_split = NULL,
  cluster_row_slices =F,
  
  top_annotation = NULL,
  column_split =NULL,
  cluster_column_slices = F,
  text_font_size = 10,
  ...){
  message('IMPORTANT: do not change the order `by` and `celltype` of res generate by `stat_metadata_counts`')
  library(ComplexHeatmap)
  
  # 1) 
  mat <- res$mat
  ComplexHeatmap::Heatmap(
    mat,
    name = 'Percentages',
    col =col,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = text_font_size))
    },
    
    # row
    right_annotation = right_annotation ,
    row_split = if(!is.null(row_split)) res$anno_by[rownames(mat),row_split[1]],
    column_split = if(!is.null(column_split)) res$anno_celltype[colnames(mat),column_split[1]] ,
    cluster_row_slices = cluster_row_slices ,
    # column
    top_annotation = top_annotation,
    cluster_column_slices = cluster_column_slices,
    ...
  )
})







