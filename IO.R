# TODO, load workspace
# wirte base load data and save data
# load/save using base load/save method for: data, analysis, dl
# write methods for params
ogel$set('public','load_data',function(load_analysis=F,data_use = 'data'){
  stopifnot(is.null(self$get_data(data_use)))
  
  cat("loading data from: ",self$params$qs,'\n')
  dta <-  qs::qread(self$params$qs,nthreads = self$threads)
  self$set_data(dta,data_use = data_use)
  if(load_analysis){
    cat("loading analysis list from: ",self$params$analysis,'\n')
    if(!file.exists(self$params$analysis)){
      message("analysis.qs not found, skip loading analysis list")
    }else{
      self$analysis <- qs::qread(self$params$analysis,nthreads = self$threads)
    }
  }
  rm(dta);gc()
  invisible(self)
})

ogel$set('public','save_data',function(save_anlaysis=F,force=F,data_use = 'data'){
  stopifnot(!is.null(self$get_data(data_use)))
  cat("saving data using : ",data_use,' slot into: ',self$params$qs,'\n')
  if(force || !file.exists(self$params$qs)){
    self$get_data(data_use) %>% qs::qsave(self$params$qs,nthreads = self$threads)
  }
  if(save_anlaysis){
    if(length(self$analysis) > 0){
      cat("saving analysis list into: ",self$params$analysis,'\n')
      qs::qsave(self$analysis,self$params$analysis,nthreads = self$threads)
    }else{
      message("analysis list is empty, skip saving analysis list")
    }
  }
  invisible(self)
})

ogel$set('public','save_analysis',function(force=F){
  if(!file.exists(self$params$analysis) || force){
    if(length(self$analysis) > 0 && !is.null(self$params$analysis)){
      cat("saving analysis list into: ",self$params$analysis,'\n')
      self$analysis %>%qs::qsave(self$params$analysis,nthreads = self$threads)
    }else{
      message("analysis list is empty, skip saving analysis list")
    }
  }
  invisible(self)
})

ogel$set('public','load_analysis',function(){
  if(file.exists(self$params$analysis)){
    if(length(self$analysis) == 0){
      self$analysis <- qs::qread(self$params$analysis,nthreads = self$threads)
    }else{
      message("analysis list is not empty, skip loading analysis list")
    }
  }else{
    message("analysis.qs not found, skip loading analysis list")
  }
  invisible(self)
})



# data
ogel$set("public","load_rds_qs",function(file_name){
  if(stringr::str_detect(file_name,stringr::regex("\\.rds$",ignore_case=TRUE))){
    invisible(base::readRDS(file_name))
  }else if(stringr::str_detect(file_name,stringr::regex("\\.qs$",ignore_case=TRUE))){
    invisible(qs::qread(file_name))
  }else{
    stop(paste0("file: ",file_name," is not a valid rds or qs file"))
  }
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



