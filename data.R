sc$set('public', 'subset', function(values = NULL, key = self$celltype,create_dir = F,drop_metadata_levels=T,tag='tmp') {

  stopifnot(!is.null(tag))
  if (is.null(values)) {
    message('No values provided, not subsetting')
    stop("No values provided, not subsetting")
  }
  idx <- self$data@meta.data[[key]] %in% values
  stopifnot(sum(idx) > 0)
  result <- self$clone(deep=T)
  message('Cells filtered: ', sum(idx))
  result$data <- result$data[, idx]
  result$tag <- tag
  result$analysis <- list()
  all_meta_features <- colnames(result$data@meta.data)
  if(drop_metadata_levels){
    cat("dropping metadata levels\n")
    for(each_feature in all_meta_features){
      if(is.factor(result$data@meta.data[[each_feature]])){
        result$data@meta.data[[each_feature]] <- droplevels(result$data@meta.data[[each_feature]])
      }
    }
  }
  if(create_dir && !dir.exists(self$path)){
    message("creating folder: ",self$path)
    dir.create(self$path,recursive = T)
  }
  invisible(result)
})

sc$set('public','subset_inplace',function(data_name,values = NULL,key = self$celltype,empty_analysis = T,tag = NULL){
  stopifnot(data_name != 'data')
  idx <- self$data@meta.data[[key]] %in% values
  stopifnot(sum(idx) > 0)
  self$data <- self$data[, idx]
  if(empty_analysis){
    self$analysis <- list()
  }
  if(!is.null(tag)){
    self$tag <- tag
  }
  invisible(self)
})

sc$set('public','add_moduleScore',function(sig_list,...,data_use = 'data'){
  for(i in seq_along(sig_list)){
    current_names <- names(sig_list)[i]
    message(current_names)
    new_name <- paste0(current_names,'1')
    
    stopifnot(!new_name %in% colnames(self$data@meta.data))
    stopifnot(!current_names %in% colnames(self$data@meta.data))
    self$data <- AddModuleScore(self$data,features =sig_list[current_names],name =current_names)
    colnames(self$data@meta.data) <-
      gsub(x = colnames(self$data@meta.data)
           , pattern = paste0(current_names,'1')
           , replacement = current_names
      )
  }
  invisible(self)
})

# get feature expression data
sc$set('public','get_features_expr_df',function(features,group_by = self$celltype,features_df=NULL,group_df=NULL,data_use = 'data'){
  stopifnot(is.null(features_df) || all(features %in% rownames(features_df))) # make sure feature
  p <- DotPlot(self[[data_use]],features = features,group.by = group_by)
  p_data <- 
    p$data %>% 
    mutate(features = features.plot,
           group_by = id,
           pct_exp = pct.exp,
           avg_expr_scaled = avg.exp.scaled
    ) %>% 
    select(features,group_by,pct_exp,avg_expr_scaled)
  # check group by
  all_groups <- unique(p_data$group_by)
  stopifnot(is.null(group_df) || all(all_groups %in% rownames(group_df))) # make sure group_by values consistent
  
  if(!is.null(features_df)){
    features_df_col <- colnames(features_df) %>% setdiff(colnames(p_data))
    features_df <- features_df[,features_df_col,drop=F]
    features_df$features <- rownames(features_df)
    p_data <- merge(p_data,features_df)
  }
  
  if(!is.null(group_df)){
    group_df_col <- colnames(group_df) %>% setdiff(colnames(p_data))
    group_df <- group_df[,group_df_col,drop=F]
    group_df$group_by <- rownames(group_df)
    p_data <- merge(p_data,group_df)
  }
  # merge data
  p_data
})


sc$set('public','get_expr2meta',function(features,assay = "RNA",inplace=F,data_use = 'data'){ 
  obj <- self$get_data(data_use)
  assay <- Seurat::GetAssayData(obj,assay = assay)[features,,drop=F]
  if(inplace){
    for(each_feature in features){
      obj@meta.data[[each_feature]] <- assay[each_feature,]
    }
    self$set_data(obj,data_use)
    invisible(self)
  }else{
    meta.data <- obj@meta.data
    for(each_feature in features){
      meta.data[[each_feature]] <- assay[each_feature,]
    }
    invisible(meta.data)
  }

})

sc$set('public','get_uniq',function(feature,data_use = 'data'){
  if(feature %in% colnames(self$data@meta.data)){
    stopifnot(!is.numeric(self$data@meta.data[[feature]]))
    return(sort(unique(self$data@meta.data[[feature]])))
  }else{
    stop(paste0("feature: ",feature," not found in meta.data"))
  } 
})



