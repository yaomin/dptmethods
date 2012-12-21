load.mappedData <-
function(smpStr, from=c("R","file"), pos=-1, pos.out=parent.frame()) {
  ##browser()
  from <- match.arg(from)
  if(from == "R") infiles <- smpStr
  else if(from == "file") infiles <- scan(get.str.smpsheet(smpStr), what="list", comment.char="#")
  files <- paste("..",infiles,sep="/")
  dataload <- list()

  for (i in seq(along=files)) {
    dataload <- c(dataload, list(read.table(files[i])))
  }
  names(dataload) <- paste("s", seq(along=dataload), sep="")
  n.reads <- unlist(lapply(dataload, function(x) sum(x$V4)))
  save(n.reads, file=file.path(get.wsoption.path("joinSample", pos=pos),"nreads.RData"))

  count.tbs <- sapply(dataload, function(x) {idf <- as.data.frame(table(x$V4));
                                             names(idf) <- c("count","freq");
                                             idf},
                      simplify=F)
  save(count.tbs, file=file.path(get.wsoption.path("joinSample", pos=pos),"countdistr.RData"))
  
  assign("dataload", dataload, pos=pos.out)
  assign("n.reads", n.reads, pos=pos.out)
}
