# write based method here, forexample, PCA, harmony, umap etc. the method in this script should be genral. all the other functionality will call this method
# PCA
ogel$set('public','pca',function(
    assay_use = 'RNA',
    dims = 1:50,
    data_use = 'data'
){
  obj <- self$get_data(data_use)
  if(inherits(obj,'Seurat')){
    cat("*   Seurat: Running PCA\n")
    obj <- RunPCA(obj,dims = dims,assay = assay_use)
  }else{
    cat("I am doing nothing, current support pca for seruat")
  }
  self$set_data(obj,data_use)
  rm(obj);gc()
  invisible(self)
})

# Preprocess
ogel$set('public','split_layers',function(assay_use = "RNA",batch=NULL,data_use = 'data'){
  cat("*   Please notice that: split layers only works for seurat object\n")
  # data prepare
  obj <- self$get_data(data_use)
  verbose <- self$verbose
  # split layers is only used for seurat
  if(inherits(obj,'Seurat')){
    if(is.null(batch)){
      cat("*   Batch not supplied, not splitting layers\n")
    }else{
      batch <- batch[1]
      batch_levels <- unique(obj@meta.data[[batch]])
      if(length(batch_levels) == 1){
        cat("*   Only one batch levels found, not splitting layers:",batch_levels,"\n")
      }else{
        cat("*   Batch levels found, splitting: ",paste(batch_levels,collapse = ","),"\n")
        obj[[assay_use]] <- split(obj[[assay_use]],obj@meta.data[[batch]])
      }
    }
  }
  self$set_data(obj,data_use)
  rm(obj);gc()
  invisible(self)
})

ogel$set('public','join_layers',function(assay_use = "RNA",data_use = 'data'){
  cat("*   Please notice that: join layers only works for seurat object\n")
  # data prepare
  obj <- self$get_data(data_use)
  verbose <- self$verbose
  # split layers is only used for seurat
  if(inherits(obj,'Seurat')){
    tryCatch({
      obj <- JoinLayers(obj,assay = assay_use)
      self$set_data(obj,data_use)
      rm(obj);gc()
    },error = function(e){
      cat("*   The layers might already be joined\n")
    })
  }
  invisible(self)
})


ogel$set('public','LogNormalize',function(
    assay_use = 'RNA',
    scale_factor = 10000,
    selection_method = "vst",
    nFeatures = 3000,
    vars_to_regress = c('percent.mt'),
    force = FALSE,
    data_use = 'data'){
  # data prepare
  obj <- self$get_data(data_use)
  verbose <- self$verbose
  if(inherits(obj,'Seurat')){
    assay_layers_found <- SeuratObject::Layers(obj,assay = assay_use)
    # "data"       "counts"     "scale.data"
    if(!'data' %in% assay_layers_found || force){
      cat("*   Seurat: Running LogNormalize\n")
      obj <- NormalizeData(obj, assay = assay_use,
                           normalization.method = "LogNormalize", 
                           scale.factor = scale_factor,
                           verbose = verbose)
    }
    if(length(SeuratObject::VariableFeatures(obj,assay = assay_use)) == 0 || force){
      cat("*   Seurat: Running FindVariableFeatures\n")
      obj <- FindVariableFeatures(obj, assay = assay_use, 
                                  selection.method = selection_method, 
                                  nfeatures = nFeatures,
                                  verbose = verbose)
    }
    if(!'scale.data' %in% assay_layers_found || force){
      cat("*   Seurat: Running ScaleData\n")
      obj <- ScaleData(obj, assay = assay_use,
                       vars.to.regress = vars_to_regress,
                       verbose = verbose)
    }
  }
  
  # set data
  self$set_data(obj,data_use)
  rm(obj);gc()
  invisible(self)
})

# SCT
ogel$set('public','SCT',function(
    assay_use = 'RNA',
    regress_vars = c('percent.mt'),
    conserve_memory = TRUE,
    force = FALSE,
    data_use = 'data'){
  stop("Not finished")
  obj <- self$get_data(data_use)
  if(inherits(obj,'Seurat')){
    cat("*   Seurat: Running SCT\n")
    obj <- SCTransform(obj, 
                       vars.to.regress = regress_vars, 
                       conserve.memory = conserve_memory)
  }
  self$set_data(obj,data_use)
  rm(obj);gc()
  invisible(self)
})

# harmony

ogel$set('public','harmony',function(
    assay_use = 'RNA',
    batch = 'Sample',
    dim_use = 1:50,
    data_use = 'data'){
    obj <- self$get_data(data_use)
    verbose <- self$verbose
    threads <- self$threads
    if(inherits(obj,'Seurat')){
        cat("*   Seurat: Running Harmony\n")
           obj <-  harmony::RunHarmony(
                    obj,
                    assay.use=assay_use, 
                    group.by.vars = batch,
                    ncores = threads,
                    verbose = verbose,
                    dims.use=dim_use)
    }
    self$set_data(obj,data_use)
    rm(obj);gc()
    invisible(self)
})


ogel$set('public','umap',function(
    reduction_use = 'harmony',
    dims = 1:30,
    reduction_name = 'umap',
    n.neighbors = 15,
    min.dist = 0.1,
    metric = 'euclidean',
    spread = 1,
    data_use = 'data'){
    obj <- self$get_data(data_use)
    verbose <- self$verbose
    threads <- self$threads

    if(inherits(obj,'Seurat')){
        cat("*   Seurat: Running UMAP\n")
           obj <-  RunUMAP(
            obj,
            reduction = reduction_use,
            reduction.name = reduction_name,
            dims =dims,
            n.neighbors = n.neighbors,
            min.dist = min.dist,
            metric = metric,
            spread = spread,
            verbose = verbose)
    }
    self$set_data(obj,data_use)
    rm(obj);gc()
    invisible(self)
})


ogel$set('public', 'find_neighbors', function(
    reduction = 'harmony',
    dims = 1:30,
    k.param = 20,
    nn.method = 'annoy',
    annoy.metric = 'euclidean',
    n.trees = 50,
    prune.SNN = 1/15,
    nn.eps = 0,
    l2.norm = TRUE,
    data_use = 'data') {
    
    obj <- self$get_data(data_use)
    verbose <- self$verbose
    threads <- self$threads

    if(inherits(obj, 'Seurat')) {
        cat("*   Seurat: Finding neighbors\n")
        obj <- FindNeighbors(
            object = obj,
            reduction = reduction,
            dims = dims,
            k.param = k.param,
            nn.method = nn.method,
            annoy.metric = annoy.metric,
            n.trees = n.trees,
            prune.SNN = prune.SNN,
            nn.eps = nn.eps,
            l2.norm = l2.norm,
            verbose = verbose
        )
    }
    self$set_data(obj, data_use)
    rm(obj); gc()
    invisible(self)
})

ogel$set('public', 'find_clusters', function(
    resolution = 0.5,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    group.singletons = TRUE,
    cluster_name = 'clusters',
    data_use = 'data') {
    
    obj <- self$get_data(data_use)
    verbose <- self$verbose

    if(inherits(obj, 'Seurat')) {
        cat("*   Seurat: Finding clusters\n")
        obj <- FindClusters(
            object = obj,
            resolution = resolution,
            algorithm = algorithm,
            n.start = n.start,
            n.iter = n.iter,
            group.singletons = group.singletons,
            cluster.name = cluster_name,
            verbose = verbose
        )
        all_cls_labels <- as.integer(obj@meta.data[[cluster_name]])
        cls_levels <- paste0('C',sort(unique(all_cls_labels)))
        obj@meta.data[[cluster_name]] <- factor(paste0("C",all_cls_labels),levels = cls_levels)
    }
    self$set_data(obj, data_use)
    rm(obj); gc()
    invisible(self)
})
