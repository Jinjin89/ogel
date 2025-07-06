# TODO, load workspace
# wirte base load data and save data
# load/save using base load/save method for: data, analysis, dl
# write methods for params
sc$set('public','load_data',function(qds_name = NULL,load_analysis=T,data_use = 'data'){
  stopifnot(is.null(self[[qds_name]]))
  stopifnot(is.null(self$get_data(data_use)))
  file_path <- file.path(self$path,qds_name)
  stopifnot(file.exists(file_path))
  
  cat("loading data from: ",file_path,'\n')
  dta <-  qs::qread(file_path,nthreads = self$threads)
  self$set_data(dta,data_use = data_use)
  if(load_analysis){
    cat("loading analysis list from: ",file.path(self$path,'analysis.qds'),'\n')
    if(!file.exists(file.path(self$path,'analysis.qds'))){
      message("analysis.qds not found, skip loading analysis list")
    }else{
      self$analysis <- qs::qread(file.path(self$path,'analysis.qds'),nthreads = self$threads)
    }
  }
  invisible(self)
})

sc$set('public','save_data',function(qds_name,save_anlaysis=T,force=F,data_use = 'data'){
  file_path <- file.path(self$path,qds_name)
  stopifnot(!is.null(self$get_data(data_use)))
  cat("save data into: ",file_path,'\n')
  cat("saving data using : ",data_use,' slot into: ',file_path,'\n')
  if(force || !file.exists(file_path)){
    self$get_data(data_use) %>% qs::qsave(file_path,nthreads = self$threads)
  }
  if(save_anlaysis){
    if(length(self$analysis) > 0){
      cat("saving analysis list into: ",file.path(self$path,paste0(qds_name,'_analysis.qs')),'\n')
      qs::qsave(self$analysis,file.path(self$path,'analysis.qds'),nthreads = self$threads)
    }else{
      message("analysis list is empty, skip saving analysis list")
    }
  }
  invisible(self)
})

sc$set('public','save_analysis',function(){
  if(length(self$analysis) > 0 && !is.null(self$params$analysis)){
    if(file.exists(self$params$analysis)){
      message("analysis.qds found, skip saving analysis list")
    }else{
      cat("saving analysis list into: ",self$params$analysis,'\n')
      self$analysis %>%qs::qsave(self$params$analysis,nthreads = self$threads)
    }
  }else{
    message("analysis list is empty, skip saving analysis list")
  }
  invisible(self)
})

sc$set('public','load_analysis',function(){
  if(file.exists(self$params$analysis)){
    if(length(self$analysis) == 0){
      self$analysis <- qs::qread(self$params$analysis,nthreads = self$threads)
    }else{
      message("analysis list is not empty, skip loading analysis list")
    }
  }else{
    message("analysis.qds not found, skip loading analysis list")
  }
  invisible(self)
})



# data
sc$set("public","load_rds_qds",function(file_name){
  if(stringr::str_detect(file_name,stringr::regex("\\.rds$",ignore_case=TRUE))){
    invisible(readRDS(file_name))
  }else if(stringr::str_detect(file_name,stringr::regex("\\.qds$",ignore_case=TRUE))){
    invisible(qs::qread(file_name))
  }else{
    stop(paste0("file: ",file_name," is not a valid rds or qds file"))
  }
})

sc$set("public","load_dl",function(file_name,dl_name=NULL,path = self$path){
  if(is.null(dl_name)) dl_name <- stringr::str_remove(file_name,"\\.rds|\\.qds$")
  cat("dl_name: ",dl_name,"\n")
  if(!is.null(self$dl[[dl_name]])){
    message("dl already loaded,skipping")
  }else{
    if(file.exists(file.path(path,file_name))){
      self$dl[[dl_name]] <- self$load_rds_qds(file.path(path,file_name))
      cat('data type: ',class(self$dl[[dl_name]]),'\n')
    }else{
      stop(paste0("file: ",file.path(path,file_name)," not found"))
    }
  }
  invisible(self)

})



