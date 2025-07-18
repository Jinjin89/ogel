ogel$set('public','log',function(...){
    cat(format(as.POSIXct(Sys.time()), format='%Y-%m-%d %H:%M:%OS '), ...)
    invisible(self)
})

ogel$set('public','sep_line',function(){
  cat('----------------------------------------\n')
})
