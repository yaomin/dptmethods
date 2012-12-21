para.testRatios <-
function(dat, x.n, y.n, w.n=2500000, cols=2:3,
                            chrs=NULL, n.clusters=2,
                            sourcefile=NULL) {
  ## dat is assumed from para.joinSamples
  dat <- lapply(dat, consolidate.lows)
  if(is.null(chrs)) chrs <- seq(length(dat))
  else if(is.character(chrs)) chrs <- match(chrs, names(dat))
  sfInit(parallel = T, cpus = n.clusters)
  sfExport('dat', local=T)
  ##sfSource(sourcefile)
  res <- sfLapply(chrs, test.ratios.chr, dat=dat, x.n=x.n, y.n=y.n, w.n=w.n, cols=cols)
  sfStop()
  names(res) <- chrs
  res
}
