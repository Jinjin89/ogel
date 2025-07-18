ogel$set('public', 'RNA_pipeline', function(
    # Data parameters
    data_use = 'data',
    assay_use = 'RNA',
    batch = 'Sample',
    
    # Preprocess
    # LogNormalize parameters
    preprocess_method = 'LogNormalize',
    force_preprocess = FALSE,

    scale_factor = 10000,
    selection_method = "vst",
    nFeatures = 3000,
    vars_to_regress = c('percent.mt'),
    
    
    # PCA parameters
    pca_dims = 1:50,
    
    # Batch correction
    # Harmony parameters
    harmony_dims = 1:50,
    
    # FindNeighbors parameters
    neighbors_reduction = 'harmony',
    neighbors_dims = 1:30,
    k_param = 20,
    nn_method = 'annoy',
    annoy_metric = 'euclidean',
    n_trees = 50,
    prune_SNN = 1/15,
    nn_eps = 0,
    l2_norm = TRUE,
    
    # FindClusters parameters
    resolution = 0.5,
    algorithm = 1,
    n_start = 10,
    n_iter = 10,
    group_singletons = TRUE,
    cluster_name = 'clusters',
    
    # UMAP parameters
    umap_reduction = 'harmony',
    umap_dims = 1:30,
    n_neighbors = 15,
    min_dist = 0.1,
    metric = 'euclidean',
    spread = 1
) {
    is_batch <- \(obj,batch){
        if(is.null(batch)){
            return(FALSE)
        }
        if(batch %in% names(obj@meta.data)){
            if(length(unique(obj@meta.data[[batch]])) > 1){
                return(TRUE)
            }
        }
        return(FALSE)
    }
    batch_index <- is_batch(self$get_data(data_use),batch)
    cat("* batch_index: ",batch_index,"\n")
    # Preprocess
    if (preprocess_method == 'LogNormalize'){
        cat("*  seurat preprocess: LogNormalize\n")
        if(batch_index && force_preprocess){
            cat("*  Is batch and lognormalize, split layers\n")
            self$split_layers(assay_use = assay_use,batch = batch,data_use = data_use)
        }
        self$LogNormalize(
            assay_use = assay_use,
            scale_factor = scale_factor,
            selection_method = selection_method,
            nFeatures = nFeatures,
            vars_to_regress = vars_to_regress,
            force = force_preprocess,
            data_use = data_use
        )
        if(batch_index && force_preprocess){
            cat("*  after lognormalization, join layers\n")
            self$join_layers(assay_use = assay_use,data_use = data_use)
        }
    } else if (preprocess_method == 'SCT') {
        stop("CURRENTLY NOT DEVELOPED")
        self$SCT(
            assay_use = assay_use,
            regress_vars = vars_to_regress,
            force = force_preprocess,
            data_use = data_use
        )
    }

    # Dimensionality reduction
    self$pca(
        dims = pca_dims,
        assay_use = assay_use,
        data_use = data_use
    )

    # Batch correction
    if(batch_index){
        cat("*  seurat batch correction: Harmony\n")
        self$harmony(
            assay_use = assay_use,
            batch = batch,
            dim_use = harmony_dims,
            data_use = data_use
        )
    }
    # Find neighbors
    reductions_found <- names(self$get_data(data_use)@reductions)
    cat("*  Find neighbors\n")
    if(!neighbors_reduction %in% reductions_found){
        if('harmony' %in% reductions_found){
            neighbors_reduction <- 'harmony'
        }else if('pca' %in% reductions_found){
            neighbors_reduction <- 'pca'
        }else{
            stop("could not find neighbors_reduction, found: ",paste0(reductions_found,collapse = ","),'supplied is: ',neighbors_reduction)
        }
        cat('*  supplied neighbors_reduction not found, use: ',neighbors_reduction,'\n')
    }else{
        cat("*  neighbors_reduction found: ",neighbors_reduction,'\n')
    }
    
    self$find_neighbors(
        reduction = neighbors_reduction,
        dims = neighbors_dims,
        k.param = k_param,
        nn.method = nn_method,
        annoy.metric = annoy_metric,
        n.trees = n_trees,
        prune.SNN = prune_SNN,
        nn.eps = nn_eps,
        l2.norm = l2_norm,
        data_use = data_use
    )

    # Clustering
    cat("*  seurat clustering\n")
    self$find_clusters(
        resolution = resolution,
        algorithm = algorithm,
        n.start = n_start,
        n.iter = n_iter,
        group.singletons = group_singletons,
        cluster_name = cluster_name,
        data_use = data_use
    )

    # UMAP
    cat("*  seurat umap\n")
    if(!umap_reduction %in% reductions_found){
        if('harmony' %in% reductions_found){
            umap_reduction <- 'harmony'
        }else if('pca' %in% reductions_found){
            umap_reduction <- 'pca'
        }else{
            stop("could not find umap_reduction, found: ",paste0(reductions_found,collapse = ","),'supplied is: ',umap_reduction)
        }
        cat('*  supplied umap_reduction not found, use: ',umap_reduction,'\n')
    }
    self$umap(
        reduction_use = umap_reduction,
        dims = umap_dims,
        n.neighbors = n_neighbors,
        min.dist = min_dist,
        metric = metric,
        spread = spread,
        data_use = data_use
    )

    invisible(self)
})
