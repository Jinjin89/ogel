library(R6)

ogel <- R6Class(
  'ogel',
  active = list(
    path = function(value){
      if(missing(value)){
        path <- file.path(self$workspace,self$tag)
        if(!dir.exists(path)) dir.create(path,recursive = T)
        path
      }else{
        stop('Can not set path, it is build using self$workspace and self$tag')
      }
    },
    meta.data = function(value){
      if(missing(value)){
        self$get_metadata()
      }else{
        stop('Can not set meta.data, it is a read only field')
      }
    }
  ),
  public = list(
    # Fields
    workspace = NULL,
    tag = NULL,
    data = NULL,
    celltype = NULL,
    group = NULL,
    threads = NULL,
    verbose = TRUE,
    analysis = list(),
    dl = list(),
    params = list(),
    db = list(),
    funs = list(),
    
    # Initialize method with celltype and group arguments
    initialize = function(workspace, tag,data=NULL,celltype = NULL, group = NULL, threads = 16,verbose = TRUE) {
      suppressMessages({
        library(future)
        library(future.apply)
        library(ggpubr)
        library(ComplexHeatmap)
        library(scplotter)
        library(qs)
        library(Seurat)
        library(SeuratObject)
        library(magrittr)
        library(dplyr)
        library(data.table)
      })
      self$workspace <- workspace
      self$tag <- tag
      self$data <- data
      self$celltype <- celltype
      self$group <- group
      self$threads <- threads
      self$verbose <- verbose
      
      # params
      self$params <- list(
        data_name = 'data.qs',
        analysis_name = 'analysis.qs',
        res_folder <- NULL,
        print_size = TRUE
      )
      if(!dir.exists(self$path)) {
        message("creating folder: ",self$path)
        dir.create(self$path, recursive = TRUE)
      }
      analysis_folder <- file.path(self$path,'analysis')
      if(!dir.exists(analysis_folder)){
        message("creating folder: ",analysis_folder)
        dir.create(analysis_folder, recursive = TRUE)
      }
      if(!is.null(self$res_folder)){
        if(!dir.exists(self$res_folder)){
          message("creating folder: ",self$res_folder)
          dir.create(self$res_folder, recursive = TRUE)
        }
      }
      invisible(self)
    },
    
    # Utility methods
    print = function() {
      cat("ogel object:\n")
      if (!is.null(self$data)) {
        cat("  path:",self$path,"\n")
        cat("  tag:",self$tag,"\n")
        cell_type <- ifelse(is.null(self$celltype),"Not supplied",self$celltype)
        group <- ifelse(is.null(self$group),"Not supplied",self$group)
        cat("  celltype:",cell_type,"\n")
        cat("  group:",group,"\n")
        cat("  Assays:", paste(names(self$data@assays), collapse = ", "), "\n")
        cat("  Active assay:", self$data@active.assay, "\n")
        cat("  Reductions:", paste(names(self$data@reductions), collapse = ", "), "\n")
        cat("  Number of cells:", ncol(self$data), "\n")
        cat("  Number of genes:", nrow(self$data), "\n")
        if(!is.null(self$celltype)){
          cat("  Cell type:", paste0(self$get_uniq(self$celltype),collapse = ","), "\n")
        }
        if(!is.null(self$group)){
          cat("  Group:", paste0(self$get_uniq(self$group),collapse = ","), "\n")
        }
      }
      if(self$params$print_size){
        cat("  data size:",self$get_data_size(),"\n")
        cat("  analysis size:",self$get_analysis_size(),"\n")
      }
      invisible(self)
    }
  )
)


# general function
ogel$set('public','get_metadata',function(data_use = 'data'){
  return(self$get_data(data_use)@meta.data)
})



ogel$set('public','set_data',function(data,data_use = 'data',data_new =NULL){
  stopifnot(data_use == 'data')
  self$data <- data
  invisible(data)
})

ogel$set('public','get_data',function(data_use = 'data'){
  stopifnot(data_use == 'data')
  self[[data_use]]
})


ogel$set('public','set_parallel',function(n_cores = self$threads,size = 500){
  plan(sequential)
  if(Sys.getenv('RSTUDIO') == '1'){
    cat("using multisession\n")
    plan(multisession,workers = self$threads)
  }else{
    cat("using multicore\n")
    plan(multicore,workers = self$threads)
  }
  cat('also set future.globals.maxSize to ',size,'GB\n')
  options(future.globals.maxSize= size*1024^3)
  invisible(self)
})
ogel$set('public','set_sequential',function(){
  plan(sequential)
  invisible(self)
})

# size info 
ogel$set("public","get_size",function(name = NULL){
  if(is.null(name)){
    format(object.size(self),units = 'auto')
  }else{
    format(object.size(self[[name]]),units = 'auto')
  }
})

ogel$set("public","get_analysis_size",function(){
  self$get_size('analysis')
})

ogel$set("public","get_data_size",function(){
  self$get_size('data')
})


len <- function(...){
  length(...)
}

sifn <- \(logical,...){
  if(!logical){
    stop(paste(...))
  }
  invisible(T)
}

.load_data_from_file <- \(file_name,nthreads = 1,...){
  if(stringr::str_detect(file_name,stringr::regex("\\.rds$",ignore_case=TRUE))){
    invisible(base::readRDS(file_name))
  }else if(stringr::str_detect(file_name,stringr::regex("\\.qs$",ignore_case=TRUE))){
    invisible(qs::qread(file_name))
  }else{
    stop(paste0("file: ",file_name," is not a valid rds or qs file"))
  }
}

.save_data_to_file <- \(data,file_name,nthreads = 1,...){
  if(stringr::str_detect(file_name,stringr::regex("\\.rds$",ignore_case=TRUE))){
    base::saveRDS(data,file_name)
  }else if(stringr::str_detect(file_name,stringr::regex("\\.qs$",ignore_case=TRUE))){
    qs::qsave(data,file_name,nthreads = nthreads)
  }else{
    stop(paste0("file: ",file_name," is not a valid rds or qs file"))
  }
}
