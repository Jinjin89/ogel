sc$set('public','load_db',function(db_folder='~/script/R6/db/',data_use = 'data'){
  all_files <- list.files(db_folder,pattern = '*.rds',full.names = T)
  cat('found [',paste0(basename(all_files),collapse = ','),'] databases from ',db_folder,'\n')
  for(each_file in all_files){
    db_name <- basename(each_file)
    message('loading db: ', db_name)
    db_name <- gsub('.rds','',db_name)
    self$db[[db_name]] <- readRDS(each_file)
  }
  invisible(self)
})
