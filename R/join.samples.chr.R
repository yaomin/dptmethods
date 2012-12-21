join.samples.chr <-
function(chr, x, ..., by=100,filePrefix=NULL, pos=-1) {
  dpt.options <- get("dpt.options", pos=pos)
  if(is.null(filePrefix)) filePrefix <- paste(get.wsoption.path("joinSample", pos=pos),
                                              dpt.options[[c("output","joinSample")]],
                                              sep="/")
  arglist <- lapply(list(x,...), function(x) subset(x, subset=(V1 == chr)))
  filename <- paste(filePrefix, chr, "RData",sep=".")
  assign(paste("s",chr, sep="."), do.call(join.samples, args=c(arglist,by=by)), pos=pos)
  ##save.env <- if(pos==-1) globalenv() else as.environment(pos)
  save(list=(paste("s",chr, sep=".")), file=filename)
}
