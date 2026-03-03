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

ogel$set("public","monocle2",function(outfile,
                                   downsample = NULL,
                                   group_by = NULL,
                                   assay = NULL, slot = "counts", expressionFamily = "negbinomial.size",
                                   features = NULL, feature_type = "HVF",
                                   disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit",
                                   max_components = 2, reduction_method = "DDRTree",
                                   norm_method = "log",
                                   residualModelFormulaStr = NULL, pseudo_expr = 1,
                                   root_state = NULL, seed = 11, data_use = 'data'){
  obj <- self$get_data(data_use)
  if (!file.exists(outfile)) {
    set.seed(seed)
    if (!"package:DDRTree" %in% search()) {
      attachNamespace("DDRTree")
    }
    suppressMessages(library(monocle))
    if(!is.null(downsample) && downsample > 0){
      message("*   downsample to(each ident): ",downsample)
      obj <- subset(obj, downsample = downsample)
      print(dim(obj))
    }
    # stopifnot('group_by should have value, if root_state is not null!' = !is.null(root_state) && !is.null(group_by))
    assay <- if (is.null(assay)) DefaultAssay(obj) else assay
    expr_matrix <- as.sparse(GetAssayData(obj, assay = assay, slot = slot))
    p_data <- obj@meta.data
    f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
    pd <- new("AnnotatedDataFrame", data = p_data)
    fd <- new("AnnotatedDataFrame", data = f_data)
    cds <- monocle::newCellDataSet(expr_matrix,
                                   phenoData = pd,
                                   featureData = fd,
                                   expressionFamily = do.call(get(expressionFamily, envir = getNamespace("VGAM")), args = list())
    )
    if (any(c("negbinomial", "negbinomial.size") %in% expressionFamily)) {
      cds <- BiocGenerics::estimateSizeFactors(cds)
      cds <- BiocGenerics::estimateDispersions(cds)
    }
    if (is.null(features)) {
      if (feature_type == "HVF") {
        features <- VariableFeatures(obj, assay = assay)
        if (length(features) == 0) {
          features <- VariableFeatures(FindVariableFeatures(obj, assay = assay), assay = assay)
        }
      }
      if (feature_type == "Disp") {
        features <- subset(monocle::dispersionTable(cds), eval(rlang::parse_expr(disp_filter)))$gene_id
      }
    }
    message("features number: ", length(features))
    cds <- monocle::setOrderingFilter(cds, features)
    p <- monocle::plot_ordering_genes(cds)

    cds <- monocle::reduceDimension(
      cds = cds,
      max_components = max_components,
      reduction_method = reduction_method,
      norm_method = norm_method,
      residualModelFormulaStr = residualModelFormulaStr,
      pseudo_expr = pseudo_expr
    )
    cds <- orderCells(cds)

    embeddings <- t(cds@reducedDimS)
    colnames(embeddings) <- paste0(cds@dim_reduce_type, "_", 1:ncol(embeddings))
    obj[[cds@dim_reduce_type]] <- CreateDimReducObject(embeddings = embeddings, key = paste0(cds@dim_reduce_type, "_"), assay = assay)
    obj[["Monocle2_State"]] <- cds[["State"]]

    if (cds@dim_reduce_type == "ICA") {
      reduced_dim_coords <- as.data.frame(t(cds@reducedDimS))
    } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
      reduced_dim_coords <- as.data.frame(t(cds@reducedDimK))
    }

    edge_df <- igraph::as_data_frame(cds@minSpanningTree)
    edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
    edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
    trajectory <- geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
    p <- Seurat::DimPlot(obj, group.by = "Monocle2_State", reduction = reduction_method, label = TRUE) +
      trajectory

    print(p)
    if (is.null(root_state)) {
      root_state <- select.list(sort(unique(cds[["State"]])), title = "Select the root state to order cells:")
      if (root_state == "" || length(root_state) == 0) {
        root_state <- NULL
      }
    }else{
      root_state <-
        pData(cds) %>%
        dplyr::count(State,!!as.name(group_by)) %>%
        magrittr::set_colnames(c("State",'group_by','n')) %>%
        dplyr::filter(group_by == root_state) %>%
        dplyr::arrange(dplyr::desc(n)) %>%
        dplyr::pull(State)
     message('root state is set to: ',root_state[1])
    }
    cds <- orderCells(cds, root_state = root_state[1])
    obj[["Monocle2_State"]] <- cds[["State"]]
    obj[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
    monocle_res <- list(cds = cds, features = features, trajectory = trajectory)
    monocle_res %>% saveRDS(outfile)
  }
  message("*   loading from rds:",outfile)
  monocle_res <- readRDS(outfile)
  obj@tools$Monocle2 <- monocle_res
  self$set_data(obj, data_use)
  invisible(self)
})
