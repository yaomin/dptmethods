para.addCounts <-
function(x, ..., n.clusters=2, chrs=NULL,
                           sourcefile=NULL) {

  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  function(x) as.character(unique(x$V1)))))
    ##chrs <- unique(unlist(lapply(list(x,...),  names)))
  }
  
  sfInit(parallel = T, cpus = n.clusters)
  sfExportAll()
  ##sfSource(sourcefile)
  
  res <- sfLapply(chrs, add.counts.chr, x, ...)
  sfStop()
  names(res) <- chrs
  res
}
