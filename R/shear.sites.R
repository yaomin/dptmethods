shear.sites <-
function(sites, gap, reads.poi, wins.sel, win.size=100) {
  stopifnot(gap>win.size)
  gap <- gap/win.size
  for (i in unique(sites$pattern)) {
    sites.i <- subset(sites, subset=pattern==i)
    ids.i <- unique(sites.i$ID)
    ## KEEP IT HERE
    ##zwins <- c(wins.sel[[as.character(i)]][,'wins'], wins.sel[['P000']][,'wins'])
    ##zwins <- c(as.integer(sites.i[,'Win']), wins.sel[['P000']][,'wins'])
    ##zwins <- as.integer(sites.i[,'Win'])
    ## KEEP IT HERE
    for ( j in ids.i) {
      sites.i.j <- subset(sites.i, subset=ID==j)
      ids.tb <- table(sites.i.j$ID)
      ids <- names(ids.tb[ids.tb>1])
      sites.ijo <- sites.i.j[order(as.numeric(sites.i.j$Win)),]
      sites.ijo.wins <- as.integer(sites.ijo$Win)
      low.ijod <- sites.ijo.wins[-1]
      up.ijod <- sites.ijo.wins[-nrow(sites.ijo)]
      sites.ijod <- data.frame(up=up.ijod,
                               low=low.ijod,
                               d=((low.ijod-up.ijod)/win.size)> gap)
      sites.ijod.sel <- subset(sites.ijod, subset=d)
      ## KEEP IT HERE before confirmation
      ## gap.TF <- apply(sites.ijod.sel,
      ##                 1,
      ##                 function(x, zw) !all(get.reads.poi(reads.poi, x[1], x[2])%in% zw), zw=zwins)
      ## KEEP IT HERE
      gap.TF <- rep(TRUE,nrow(sites.ijod.sel))
      
      sites.ijod.sel.gap <- cbind(sites.ijod.sel, gap=gap.TF)
      if(nrow(sites.ijod.sel.gap) > 0) {
        ids.ijo <- sites.ijo$ID
        ids.ijo.out <- ids.ijo
        wins.ijo <- as.integer(sites.ijo$Win)
        aid <- 1
        
        for ( k in seq(nrow(sites.ijod.sel.gap))) {
          
          gap.i <- sites.ijod.sel.gap[k,'gap']
          gap.w <- sites.ijod.sel.gap[k,'up']
          if(gap.i) {
            ids.ijo.out[wins.ijo > gap.w] <- paste(ids.ijo[wins.ijo > gap.w], aid, sep=".")
            aid <- aid + 1
          }
        }

        sites.ijo$ID <- ids.ijo.out
        sites[sites$pattern==i & sites$ID ==j,] <- sites.ijo
      }
    }
  }
  sites
}
