# TODO, load workspace
# wirte base load data and save data
# load/save using base load/save method for: data, analysis, dl
# write methods for params
ogel$set('public','load_data',function(load_analysis=F,data_use = 'data'){
  data_file <- file.path(self$path,self$params$data_name)
  cat("loading data from: ",data_file,'\n')
  dta <-.load_data_from_file(file_name = data_file,nthreads = self$threads)
  self$set_data(dta,data_use = data_use)
  if(load_analysis){
    self$load_analysis()
  }
  rm(dta);gc()
  invisible(self)
})

ogel$set('public','save_data',function(save_anlaysis=F,force=F,data_use = 'data'){
  stopifnot(!is.null(self$get_data(data_use)))
  data_file <- file.path(self$path,self$params$data_name)
 if(force || !file.exists(data_file)){
    cat("saving data using : ",data_use,' slot into: ',data_file,'\n')
    self$get_data(data_use) %>% .save_data_to_file(data_file,nthreads = self$threads)
  }else{
    cat("File found, skip saving data\n")
  }
  if(save_anlaysis){
    self$save_analysis(force = force)
  }
  invisible(self)
})

ogel$set('public','save_analysis',function(force=F){
  analysis_file <- file.path(self$path,'analysis',self$params$analysis_name)
  if(force || !file.exists(analysis_file)){
    if(length(self$analysis) > 0 && !is.null(analysis_file)){
      cat("saving analysis list into: ",analysis_file,'\n')
      if(!dir.exists(file.path(self$path,'analysis'))){
        dir.create(file.path(self$path,'analysis'),recursive = T)
      }
      self$analysis %>% .save_data_to_file(analysis_file,nthreads = self$threads)
    }else{
      message("analysis list is empty, skip saving analysis list")
    }
  }else{
    cat("analysis results found, skip save data!")
  }
  invisible(self)
})

ogel$set('public','load_analysis',function(force = F){
  analysis_file <- file.path(self$path,'analysis',self$params$analysis_name)
  if(file.exists(analysis_file)){
    if(force || length(self$analysis) == 0){
      cat("loading analysis results")
      self$analysis <- .load_data_from_file(analysis_file,nthreads = self$threads)
    }else{
      message("analysis list is not empty, skip loading analysis list")
    }
  }else{
    message("analysis.qs not found, skip loading analysis list")
  }
  invisible(self)
})


ogel$set("public","load_dl",function(file_name,dl_name=NULL,path = self$path){
  if(is.null(dl_name)) dl_name <- stringr::str_remove(file_name,"\\.rds|\\.qs$")
  cat("dl_name: ",dl_name,"\n")
  if(!is.null(self$dl[[dl_name]])){
    message("dl already loaded,skipping")
  }else{
    if(file.exists(file.path(path,file_name))){
      self$dl[[dl_name]] <- self$load_rds_qs(file.path(path,file_name))
      cat('data type: ',class(self$dl[[dl_name]]),'\n')
    }else{
      stop(paste0("file: ",file.path(path,file_name)," not found"))
    }
  }
  invisible(self)
  
})



