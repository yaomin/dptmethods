search.sites <-
function(wins.sel, events.test, chr, win.size=100, den.cut=5e5){
  sites <- NULL
  for.seq <- match.mx(t(events.test), t(attr(wins.sel, 'event')))
  site.names <- names(wins.sel)
  for (i in for.seq) {
    den.i <- win.den(wins.sel[[i]],
                     win.size=win.size)
    sites.i <- win.search(den.i, cut=den.cut, win.size=win.size)
    chr.i <- rep(chr, nrow(sites.i))
    patt.i <- rep(site.names[i], nrow(sites.i))
    
    sites <- rbind(sites, cbind(chr.i, patt.i, sites.i))
  }
  names(sites) <- c("chr", "pattern", names(sites.i))
  attr(sites, 'event') <- events.test
  sites
}
