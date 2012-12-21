para.joinCounts <-
function(x, ..., n.clusters=2, chrs=NULL,
                            sourcefile=NULL,
                            winsize=100) {
  
  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  names)))
  }
  
  sfInit(parallel = T, cpus = n.clusters)
  sfExportAll()
  ##sfSource(sourcefile)
  res <- sfLapply(chrs, join.counts.chr, x, ...,by=winsize)
  sfStop()
  
  names(res) <- chrs
  res
}
