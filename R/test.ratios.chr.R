test.ratios.chr <-
function(chr, dat, x.n, y.n, w.n=2500000, cols=2:3) {
  ## test only on one chr
  
  dat.sub <- dat[[chr]][,cols]
  res <- as.data.frame(test.ratios(dat.sub, x.n, y.n, w.n))
  names(res) <- c('pval','ratio','rate1', 'rate2', 'ciL','ciU')
  res
}
