# multiome preprocessing
ogel$set('public','multi_preprocess',function(
    RNA_assay = 'RNA', ATAC_assay = 'peaks',
    force_rna = FALSE, force_atac = FALSE,
    regress_vars = c('percent.mt'), min_cutoff = 5, 
    conserve_memory = TRUE, data_use = 'data'){
  
  message("*   Starting multiome preprocessing")
  message("*   Processing RNA and ATAC-seq assays")

  self$multi_preprocess_rna(RNA_assay = RNA_assay, regress_vars = regress_vars, conserve_memory = conserve_memory, force = force_rna, data_use = data_use)
  self$multi_preprocess_atac(ATAC_assay = ATAC_assay, min_cutoff = min_cutoff, force = force_atac, data_use = data_use)
  message("*   Multiome preprocessing completed, now you many need to run multi_harmony")
  invisible(self)
})

# RNA preprocessing
ogel$set('public','multi_preprocess_rna',function(
    RNA_assay = 'RNA',
    regress_vars = c('percent.mt'),
    conserve_memory = TRUE,
    force = FALSE,
    data_use = 'data'){
  
  message("*   Starting RNA preprocessing")
  
  # Get the data object
  obj <- self$get_data(data_use)

  # check SCT
  if('SCT' %in% names(obj@assays) && !force) {
    message("*   SCT assay already exists, skipping SCTransform")
    invisible(self)
  }

    # Set RNA assay as default
  DefaultAssay(obj) <- RNA_assay
  
  # SCTransform with regression
  message("*   Running SCTransform...")
  obj <- SCTransform(obj, 
                    vars.to.regress = regress_vars, 
                    conserve.memory = conserve_memory)
  
  # Run PCA
  message("*   Running PCA...")
  obj <- RunPCA(obj)
  
  # Store results
  self$set_data(obj, data_use)
  message("*   RNA preprocessing completed")
  
  invisible(self)
})


# ATAC preprocessing
ogel$set('public','multi_preprocess_atac',function(
    ATAC_assay = 'peaks',
    min_cutoff = 5,
    force = FALSE,
    data_use = 'data'){

  message("*   Starting ATAC preprocessing")

  # Get the data object
  obj <- self$get_data(data_use)

  # check LSI
  if('lsi' %in% names(obj@reductions) && !force) {
    message("*   LSI reduction already exists, skipping LSI")
    invisible(self)
  }

  # ===== ATAC Processing =====
  message("*   Processing ATAC-seq assay...")
  DefaultAssay(obj) <- ATAC_assay
  
  # ATAC-seq preprocessing
  obj <- FindTopFeatures(obj, min.cutoff = min_cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  
  # Store results
  self$set_data(obj, data_use)
  message("*   ATAC preprocessing completed")

  invisible(self)

})

# Harmony batch correction
ogel$set('public','multi_harmony',function(batch_var = 'sample', RNA_assay = 'SCT', 
                                         force_lsi = FALSE,force_pca = FALSE,
                                         ATAC_assay = 'peaks', data_use = 'data'){
  
  message("*   Starting Harmony batch correction for RNA and ATAC")
  
  # Get the data object
  obj <- self$get_data(data_use)
  
  # ===== Harmony Integration =====
  
  # Harmony on LSI (ATAC)
  if('harmony_lsi' %in% names(obj@reductions) && !force_lsi) {
    message("*   Harmony LSI reduction already exists, skipping Harmony")
  }else{
    message("*   Running Harmony on ATAC (LSI) reduction...")
  obj <- RunHarmony(
    object = obj,
    group.by.vars = batch_var,
    reduction = 'lsi',
    assay.use = ATAC_assay,
    project.dim = FALSE,
    reduction.save = "harmony_lsi"
  )
  }
  
  if('harmony_pca' %in% names(obj@reductions) && !force_pca) {
    message("*   Harmony PCA reduction already exists, skipping Harmony")
  }else{
    message("*   Running Harmony on RNA (PCA) reduction...")
  obj <- RunHarmony(
    object = obj,
    group.by.vars = batch_var,
    reduction = 'pca',
    assay.use = RNA_assay,
    project.dim = FALSE,
    reduction.save = "harmony_pca"
  )
  }
  
  # Store results
  self$set_data(obj, data_use)
  message("*   Harmony batch correction completed, now you many need to run multi_integration")
  
  invisible(self)
})

ogel$set('public','multi_integration',function(rna_dims = 1:30, atac_dims = 2:30, 
                                             RNA_assay = 'SCT', reduction.name = 'umap_merge', 
                                             resolution = 0.8,
                                             data_use = 'data'){
  
  message("*   Starting multimodal integration")
  
  # Get the data object
  obj <- self$get_data(data_use)
  
  # ===== Multimodal Integration =====
  message("*   Finding multimodal neighbors...")
  DefaultAssay(obj) <- RNA_assay
  
  obj <- FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("harmony_pca", "harmony_lsi"), 
    dims.list = list(rna_dims, atac_dims),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
  )

  # ===== FindClusters =====
  message("*   Finding clusters...")
  obj <- FindClusters(
    object = obj,
    resolution = resolution
  )

  # ===== UMAP Embedding =====
  message("*   Computing UMAP embedding with name: ", reduction.name)
  obj <- RunUMAP(
    object = obj,
    nn.name = "weighted.nn",
    verbose = TRUE,
    reduction.name = reduction.name
  )
  
  # Store results
  self$set_data(obj, data_use)
  message("*   Multimodal integration completed")
  
  invisible(self)
})

ogel$set('public','multi_pipeline',function(RNA_assay = 'RNA', ATAC_assay = 'peaks',
                                          force_rna = FALSE, force_atac = FALSE,
                                          force_lsi = FALSE, force_pca = FALSE,
                                          rna_dims = 1:30, atac_dims = 2:30, 
                                          batch_var = 'sample', regress_vars = c('percent.mt'),
                                          min_cutoff = 5, conserve_memory = TRUE, 
                                          reduction.name = 'umap_merge', data_use = 'data'){
  
  message("*   Running complete multiome pipeline")
  
  # Run all three steps sequentially
  self$multi_preprocess(RNA_assay = RNA_assay, ATAC_assay = ATAC_assay,
                        force_rna = force_rna, force_atac = force_atac,
                        regress_vars = regress_vars, min_cutoff = min_cutoff,
                        conserve_memory = conserve_memory, data_use = data_use)
  
  self$multi_harmony(batch_var = batch_var, RNA_assay = 'SCT', 
                     ATAC_assay = ATAC_assay, force_lsi = force_lsi, 
                     force_pca = force_pca, data_use = data_use)
  
  self$multi_integration(rna_dims = rna_dims, atac_dims = atac_dims,
                         RNA_assay = 'SCT', reduction.name = reduction.name, 
                         data_use = data_use)
  
  message("*   Complete multiome pipeline finished successfully")
  
  invisible(self)
})

ogel$set('public','multi_macs2',function(...,atac_assay = 'peaks',cleanup=FALSE,macs2.path = '/opt/conda/envs/macs2/bin/macs3', data_use = 'data'){
  message("*   Running MACS2 peak calling")
  
  # Get the data object
  obj <- self$get_data(data_use)
  obj@active.assay <- atac_assay
  peaks_dir <- file.path(self$path, 'peaks')
  if(!dir.exists(peaks_dir)) {
    dir.create(peaks_dir,recursive = T)
  }
  # Run callPeaks
  self$analysis$peaks <- CallPeaks(obj,outdir = peaks_dir, macs2.path = macs2.path,cleanup=cleanup, ...)
  self$analysis$peaks  %>% qs::qsave(file.path(self$path, 'macs2peaks.qds'),nthreads = self$threads)
  invisible(self)
})

# TODO read peaks, or load atac data from raw data