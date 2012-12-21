map.region.2 <-
function(tst, gap=500) {
  chrs <- unique(tst$chr)
  out.1 <- NULL
  for(chrom in chrs) {
    tst.chr <- subset(tst, subset=chr==chrom)
    out.1 <- rbind(out.1, map.region.chr.2(tst.chr, gap))
  }
  
  r.range <- t(sapply(by(out.1[,c('start','end')], out.1$R.id, range), function(x) x))
  out.2 <- cbind(out.1,r.range[match(out.1$R.id, row.names(r.range)),])
  names(out.2) <- c(names(out.1), c('r.start','r.end'))
  out.2
}
