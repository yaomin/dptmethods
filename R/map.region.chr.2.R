map.region.chr.2 <-
function(tst, gap=500) {
 
  tst.o <- tst[order(as.numeric(tst$start)),]
  tst.d <- tst.o$start[-1]-tst.o$end[-nrow(tst.o)]
  tst.flw <- tst.d < gap

  tst.id <- rep(0, nrow(tst))
  for(j in seq(nrow(tst.o))) {
    if(j==1) {
      tst.id[j] <- 1
    }
    else {
      if(tst.flw[j-1]) tst.id[j] <- tst.id[j-1]
      else {
        tst.id[j] <- tst.id[j-1] + 1
      }
    }
  }

  r.id <- paste("R", tst.o$chr, tst.id, sep="-")

  transform(tst.o, R.id=r.id)
  
  
}
