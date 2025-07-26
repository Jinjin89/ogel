ogel$set("public","slingshot",function(group_by = self$celltype, reduction = 'umap', prefix = NULL, 
                                    dims = NULL, start = NULL, end = NULL, 
                                    align_start = FALSE, reverse = FALSE, seed = 11, ...,data_use = 'data'){
  suppressMessages(library(slingshot))
  if (missing(group_by)) {
    stop("group_by is missing")
  }
  if (is.null(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }
  obj <- self$get_data(data_use)
  obj_sub <- obj[, !is.na(obj[[group_by, drop = TRUE]])]
  if (is.null(dims)) {
    dims <- 1:2
  }
  
  set.seed(seed)
  sl <- slingshot(
    data = as.data.frame(obj_sub[[reduction]]@cell.embeddings[, dims]),
    clusterLabels = as.character(obj_sub[[group_by, drop = TRUE]]),
    start.clus = start, end.clus = end, ...
  )
  
  obj@tools[[paste("Slingshot", group_by, reduction, sep = "_")]] <- sl
  df <- as.data.frame(slingPseudotime(sl))
  colnames(df) <- paste0(prefix, colnames(df))
  if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
      df <- max(df, na.rm = TRUE) - df
    }
  }
  obj <- AddMetaData(obj, metadata = df)
  obj <- AddMetaData(obj, metadata = slingBranchID(sl), col.name = paste0(prefix, "BranchID"))
  self$set_data(obj,data_use)
  invisible(self)
})
