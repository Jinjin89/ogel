library(R6)

sc <- R6Class(
  'sc',
  active = list(
    path = function(value){
      if(missing(value)){
        path <- file.path(self$workspace,self$tag)
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
    res_folder = NULL,
    analysis = list(),
    dl = list(),
    params = list(),
    db = list(),
    funs = list(),

    # Initialize method with celltype and group arguments
    initialize = function(workspace, tag,data=NULL,celltype = NULL, group = NULL, threads = 16) {
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
      self$res_folder <- NULL

      # params
      self$params <- list(
        qds = file.path(self$path,'data.qds'),
        analysis = file.path(self$path,'analysis.qds'),
        analysis_folder = file.path(self$path,'analysis'),
        print_size = TRUE
      )
      if(!dir.exists(self$path)) {
        message("creating folder: ",self$path)
        dir.create(self$path, recursive = TRUE)
      }
      if(!dir.exists(self$params$analysis_folder)){
        message("creating folder: ",self$params$analysis_folder)
        dir.create(self$params$analysis_folder, recursive = TRUE)
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
      cat("sc object:\n")
      cat("  Cell type:", self$celltype %||% "Not specified", "\n")
      cat("  Group:", self$group %||% "Not specified", "\n")
      if (!is.null(self$data)) {
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

# Helper function for NULL handling
`%||%` <- function(x, y) if (is.null(x)) y else x

# general function
sc$set('public','get_metadata',function(data_use = 'data'){
  return(self$get_data(data_use)@meta.data)
})



sc$set('public','set_data',function(data,data_use = 'data',data_new =NULL){
  stopifnot(data_use == 'data')
  self$data <- data
  invisible(data)
})

sc$set('public','get_data',function(data_use = 'data'){
  stopifnot(data_use == 'data')
  self[[data_use]]
})

sc$set('public','sep_line',function(){
  cat('----------------------------------------\n')
})

sc$set('public','set_parallel',function(n_cores = self$threads,size = 500){
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
sc$set('public','set_sequential',function(){
  plan(sequential)
  invisible(self)
})

# size info 
sc$set("public","get_size",function(name = NULL){
  if(is.null(name)){
    format(object.size(self),units = 'auto')
  }else{
    format(object.size(self[[name]]),units = 'auto')
  }
})

sc$set("public","get_analysis_size",function(){
    self$get_size('analysis')
})

sc$set("public","get_data_size",function(){
    self$get_size('data')
})

# preprocess method
source('analysis.R')
source('method.R')
source('pp.R')
source('plot.R')
source('data.R')
source('multiome.R')
source('diff.R')
source('IO.R')
source('db.R')


len <- function(...){
  length(...)
}

sifn <- \(logical,...){
  if(!logical){
    stop(paste(...))
  }
  invisible(T)
}