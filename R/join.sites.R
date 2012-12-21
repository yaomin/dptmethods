join.sites <-
function(sites, gap=500, win.size=100) {
  ##TODO: a c version
  stopifnot(gap>win.size)
  gap <- gap/win.size
  patts <- unique(sites$pattern)
  res <- sites
  for (i in patts) {
    sub.i <- which(sites$pattern == i)
    sites.i <- sites[sub.i,]
    sites.by <- by(as.numeric(sites.i$Win),sites.i$ID, range)

    sites.by.i <- as.data.frame(matrix(unlist(sites.by), ncol=2, byrow=T))
    names(sites.by.i) <- c('start','end')
    row.names(sites.by.i) <- dimnames(sites.by)[[1]]

    sites.by.i.o <- sites.by.i[order(sites.by.i[,1]),]
    sites.by.i.o.IDs <- row.names(sites.by.i.o)
    sites.by.i.d <- data.frame(lw=sites.by.i.o.IDs[-1],
                               up=sites.by.i.o.IDs[-nrow(sites.by.i)],
                               d=((sites.by.i.o[-1,1]-sites.by.i.o[-nrow(sites.by.i),2])/win.size)<= gap)
    sites.by.i.d.sel <- sites.by.i.d[which(sites.by.i.d$d),]

    if (nrow(sites.by.i.d.sel) > 0) {
      tst.id <- as.integer(row.names(sites.by.i.d.sel))
      tst.lw <- as.character(sites.by.i.d.sel$lw)
      tst.up <- as.character(sites.by.i.d.sel$up)

      if(nrow(sites.by.i.d.sel) > 1) {
        for ( i in 2:length(tst.id)) {
          if((tst.id[i]-tst.id[i-1])==1) tst.up[i] <- tst.up[i-1]
        }
      }
    
      for ( j in seq(length(tst.lw))) {
        sites.i$ID[sites.i$ID==tst.lw[j]] <- tst.up[j]
      }
    }
    
    res[sub.i,] <- sites.i
  }
  res
}
