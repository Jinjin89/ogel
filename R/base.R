library(R6)

ogel <- R6Class(
  'ogel',
  private = list(
    R_DIR = this.path::this.dir(),
    DB_DIR = file.path(dirname(this.path::this.dir()),"db"),
    custom_fig_dir = NULL,
    custom_watermark_dir = NULL
  ),
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
    },
    fig_dir = function(value){
      if(missing(value)){
        # Use custom path if set, otherwise default
        if(!is.null(private$custom_fig_dir)) {
          fig_path <- private$custom_fig_dir
        } else {
          fig_path <- file.path(self$path, 'figs')
        }
        if(!dir.exists(fig_path)) dir.create(fig_path, recursive = T)
        fig_path
      }else{
        stop('Can not set fig_dir directly, use set_fig_dir() method instead')
      }
    },
    watermark_dir = function(value){
      if(missing(value)){
        # Use custom path if set, otherwise default
        if(!is.null(private$custom_watermark_dir)) {
          watermark_path <- private$custom_watermark_dir
        } else {
          watermark_path <- file.path(self$path, 'watermark')
        }
        if(!dir.exists(watermark_path)) dir.create(watermark_path, recursive = T)
        watermark_path
      }else{
        stop('Can not set watermark_dir directly, use set_watermark_dir() method instead')
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
    res_folder = NULL,
    log_file = NULL,
    analysis = list(),
    dl = list(),
    params = list(),
    plot_config = list(),
    db = list(),
    funs = list(),
    
    # Setup plot defaults - simple plot configuration
    setup_plot_defaults = function(auto_watermark = TRUE, fig_dir = NULL, watermark_dir = NULL, plot = TRUE, fontfamily = "ArialMT", fontsize = 14, ensure_dirs = FALSE) {
      # Set directories if provided
      if(!is.null(fig_dir)) {
        if(!dir.exists(fig_dir)) {
          dir.create(fig_dir, recursive = TRUE)
          message("Created fig directory: ", fig_dir)
        }
        private$custom_fig_dir <- fig_dir
      }
      if(!is.null(watermark_dir)) {
        if(!dir.exists(watermark_dir)) {
          dir.create(watermark_dir, recursive = TRUE)
          message("Created watermark directory: ", watermark_dir)
        }
        private$custom_watermark_dir <- watermark_dir
      }
      
      # Auto-create required directories when ensure_dirs is TRUE
      if(ensure_dirs) {
        # Trigger creation of fig_dir by accessing the active field
        fig_path <- self$fig_dir
        # Trigger creation of watermark_dir by accessing the active field  
        watermark_path <- self$watermark_dir
      }
      
      # Set plot config values
      self$plot_config <- list(
        auto_watermark = auto_watermark,
        plot = plot,
        fontfamily = fontfamily,
        fontsize = fontsize
      )
      
      invisible(self)
    },
    
    # Initialize method with celltype and group arguments
    initialize = function(workspace, tag,data=NULL,celltype = NULL, group = NULL, threads = 16,verbose = TRUE, log_file = NULL) {
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
        library(stringr)
      })
      self$workspace <- workspace
      self$tag <- tag
      self$data <- data
      self$celltype <- celltype
      self$group <- group
      self$threads <- threads
      self$verbose <- verbose
      
      # Handle log_file path
      if (!is.null(log_file)) {
        if (dirname(log_file) == ".") {
          # If it's just a filename, put it in workspace/tag directory
          self$log_file <- file.path(workspace, tag, log_file)
        } else {
          # If it's a path, use as is
          self$log_file <- log_file
        }
      } else {
        self$log_file <- NULL
      }
      
      
      # params
      self$params <- list(
        data_name = 'data.qs',
        analysis_name = 'analysis.qs',
        res_folder = NULL,
        print_size = TRUE
      )
      
      # Initialize plot config with defaults
      self$setup_plot_defaults(ensure_dirs = FALSE)
      
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
      # set private args
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
      cat("  R_DIR:",private$R_DIR,"\n")
      cat("  DB_DIR:",private$DB_DIR,"\n")
      
      invisible(self)
    }
  )
)


# Logging function
ogel$set('public','say',function(message, prefix = "*   "){
  full_message <- paste0(prefix, message)
  
  # Always print to console
  cat(full_message, "\n")
  
  # Also write to log file if specified
  if(!is.null(self$log_file)){
    # Add timestamp for log file
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_entry <- paste0("[", timestamp, "] ", full_message, "\n")
    
    # Write to log file
    cat(log_entry, file = self$log_file, append = TRUE)
  }
  
  invisible(self)
})

# Announcement/Header function
ogel$set('public','announce',function(message, sep_char = "=", sep_length = 60){
  # Create separator line
  separator <- paste(rep(sep_char, sep_length), collapse = "")
  
  # Get timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  time_line <- paste0("TIME: ", timestamp)
  
  # Center align the text
  center_text <- function(text, width) {
    text_length <- nchar(text)
    if (text_length >= width) return(text)
    
    padding <- width - text_length
    left_pad <- floor(padding / 2)
    right_pad <- ceiling(padding / 2)
    
    paste0(paste(rep(" ", left_pad), collapse = ""), text, paste(rep(" ", right_pad), collapse = ""))
  }
  
  # Center align both lines
  centered_message <- center_text(message, sep_length)
  centered_time <- center_text(time_line, sep_length)
  
  # Print to console
  cat("\n", separator, "\n", sep = "")
  cat(centered_message, "\n")
  cat(centered_time, "\n")
  cat(separator, "\n\n", sep = "")
  
  # Also write to log file if specified
  if(!is.null(self$log_file)){
    log_entry <- paste0("\n", separator, "\n", centered_message, "\n", centered_time, "\n", separator, "\n\n")
    cat(log_entry, file = self$log_file, append = TRUE)
  }
  
  invisible(self)
})

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
  if(is.character(data_use)){
    if(data_use == 'data'){
      return(self$data)
    }else{
      stop('data_use is not a valid data_use')
    }
  }else{
    cat('data_use is not a character, returning it directly\n')
    return(data_use)
  }

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

# Directory setting functions
ogel$set('public','set_fig_dir',function(path = NULL){
  if(is.null(path)) {
    private$custom_fig_dir <- NULL
    self$say("Fig directory reset to default")
  } else {
    if(!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      self$say(paste("Created fig directory:", path))
    }
    private$custom_fig_dir <- path
    self$say(paste("Fig directory set to:", path))
  }
  invisible(self)
})

ogel$set('public','set_watermark_dir',function(path = NULL){
  if(is.null(path)) {
    private$custom_watermark_dir <- NULL
    self$say("Watermark directory reset to default")
  } else {
    if(!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      self$say(paste("Created watermark directory:", path))
    }
    private$custom_watermark_dir <- path
    self$say(paste("Watermark directory set to:", path))
  }
  invisible(self)
})

# Plot configuration methods
ogel$set('public','set_plot_config',function(fontfamily = NULL, fontsize = NULL, auto_watermark = NULL, plot = NULL){
  if(!is.null(fontfamily)) {
    self$plot_config$fontfamily <- fontfamily
    self$say(paste("Font family set to:", fontfamily))
  }
  if(!is.null(fontsize)) {
    self$plot_config$fontsize <- fontsize
    self$say(paste("Font size set to:", fontsize))
  }
  if(!is.null(auto_watermark)) {
    self$plot_config$auto_watermark <- auto_watermark
    self$say(paste("Auto watermark set to:", auto_watermark))
  }
  if(!is.null(plot)) {
    self$plot_config$plot <- plot
    self$say(paste("Plot setting set to:", plot))
  }
  invisible(self)
})

ogel$set('public','get_plot_config',function(){
  return(self$plot_config)
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
    invisible(qs::qread(file_name,nthreads = nthreads))
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
