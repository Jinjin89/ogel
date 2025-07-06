sc$set('public','sc_pipeline',function(batch = self$sample, nfeatures=3000, data_use = 'data', data_new=data_use){
  pp_file_name <- 'pp.qds'
  file_path <- file.path(self$path,pp_file_name)
  if(!file.exists(file_path)) {
    message("*   input data should be filtered by QC!")
    message("*   input layers should be split")
    message("*   use harmony by default")
    
    # Use local variable to avoid confusion
    #obj <- do.call(fun_data,list(self))
    obj <- self$get_data(data_use)
    obj <- JoinLayers(obj)
    obj <- split(obj, f = obj[[batch]][[1]])
    # normalize
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    
    # harmony
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "harmony",
      verbose = FALSE
    )
    # clusters
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
    obj <- FindClusters(obj, cluster.name = "harmony_clusters")
    
    # umap
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap")
    
    # join layer
    obj <- JoinLayers(obj)
    
    # save the results
    message("*   saving file into: ", file_path)
    sc$set_data(obj,data_use)
    rm(obj);gc()
    self$save_data(pp_file_name,data_use = data_new)
  }
  message("*   loading file from: ", file_path)
  self$load_data(pp_file_name,data_new)
  invisible(self)
})


sc$set('public','pp_umap',function(...,data_use = 'data'){
  # TODO: improve this support modify in place
  obj <- RunUMAP(self$get_data(data_use), ...)
  self$set_data(obj,data_use)
  invisible(self)
})
