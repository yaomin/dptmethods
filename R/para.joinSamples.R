para.joinSamples <-
function(x, ..., n.clusters=2, chrs=NULL, by=100, filePrefix=NULL, pos=-1) {
  dpt.options <- get("dpt.options", pos=pos)
  if(is.null(filePrefix)) filePrefix <- paste(get.wsoption.path("joinSample", pos=pos),
                                              dpt.options[[c("output","joinSample")]],
                                              sep="/")
  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  function(x) as.character(unique(x$V1)))))
  }
  
  res <- sfLapply(chrs, join.samples.chr, x, ..., filePrefix=filePrefix,by=by)
}
