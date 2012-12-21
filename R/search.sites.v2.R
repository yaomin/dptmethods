search.sites.v2 <-
function(wins.sel,
                            events.test=NULL,
                            chr,
                            win.size=100,
                            den.cutoff = c(1e5, 5e4),
                            cut.1st=c(3000,2000))
{
### omit empty wins.sel
  noemp <- sapply(wins.sel, function(x) nrow(x)>0)
  attr.event <- attr(wins.sel, 'event')[,noemp,drop=F]
  attr.cutoff <- attr(wins.sel, 'cutoff')[noemp]
  wins.sel <- wins.sel[noemp]
  attr(wins.sel, 'event') <- attr.event
  attr(wins.sel, 'cutoff') <- attr.cutoff
 
  if(is.null(events.test)) for.seq <- seq(wins.sel)
  else {
    wins.sel <- pattern.events(wins.sel, events.test)
    for.seq <- seq(wins.sel)
  }
  
  den.cutoff <- den.cutoff[attr(wins.sel,'nonull'),,drop=F]
  events.test <- events.test[,attr(wins.sel,'nonull'),drop=F]
  ##for.seq <- match.mx(t(events.test), t(attr(wins.sel, 'event')))

  site.names <- names(wins.sel)
  if(is.null(nrow(den.cutoff))) den.cutoff <- matrix(rep(den.cutoff,
                                                         each=length(for.seq)),
                                                     ncol=2)
  else if(nrow(den.cutoff) != length(for.seq)) stop("cutoff dim does not match event dim!")

  sfExport("wins.sel","events.test","chr","win.size","den.cutoff","cut.1st")
  wins.sel.cuts <- sfSapply(for.seq,
                            function(x) cut.wins(wins.sel[[x]], den.cutoff[x,1],cut.1st[1],cut.1st[2]),
                            simplify=F)
  wins.sel.cuts.size <- sapply(wins.sel.cuts,
                               seq,
                               simplify=F)
  idx <- NULL
  for(i in seq(along=wins.sel.cuts)) {
    idx <- rbind(idx, cbind(i, seq(length(wins.sel.cuts[[i]]))))
  }
  idx.list <- as.list(as.data.frame(t(idx)))
  process.search <- function(l) {
    ##browser()
    den.cut.2 <- den.cutoff[l[1],2]
    wins.sel.cut <- wins.sel.cuts[[l[1]]][[l[2]]]
    if(nrow(wins.sel.cut)>1) {
      den.i <- win.den(wins.sel.cut,win.size=win.size)
      sites.i <- win.search(den.i, cut=den.cut.2,win.size=win.size)
    }
    else if(nrow(wins.sel.cut) == 1) {
      sites.i <- data.frame(ID=1, Win=wins.sel.cut$wins, Pr=0)
    }
    else stop('nrow(wins.sel.cut) == 0')
    sites.i$ID <- paste(l[2], sites.i$ID, sep=".")
    chr.i <- rep(chr, nrow(sites.i))
    patt.i <- rep(site.names[l[1]], nrow(sites.i))
    cbind(chr.i, patt.i, sites.i)
  }
  ##sfSource(src.file)
  
  out.list <- sfClusterApplyLB(idx.list, process.search)
  out <- NULL
  for (i in seq(along=out.list)) {
    out <- rbind(out, out.list[[i]])
  }
  names(out) <- c("chr", "pattern", names(out)[-c(1,2)])
  attr(out, 'patterns') <- events.test
  attach.winScoreP(out, wins.sel)
}
