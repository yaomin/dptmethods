#  @export

bf.filter <-
function(e.res,
                      sites.js.ex,
                      winsize,
                      bf.cutoff,
                      reads.poi,
                      pmeans.norm,
                      events.wins,
                      sample.select,
                      index.col=c("chr","pattern","ID"),
                      summary.fun=mean, sample.size=10000,
                      which.score=c("score.1","score","score.2","score.P","qvalue")) {
  ##browser()
  ## score.1 seems to be the most reasonable option
  if(missing(sample.select)) sample.select <- rep(T, ncol(pmeans.norm))
  which.score <- match.arg(which.score)
  sites.js.bfcut.1 <- mk.sites.report(sites.js.ex,
                                    index.col,
                                    "Win","Win.score",c("max","mean","min"),
                                    winsize)
  
  ##e.res.null <- e.res.null(e.res, reads.poi,sites.js.ex)

  if(all(dim(events.wins)==1)) {
    null.idx <- setdiff(row.names(pmeans.norm), sites.js.ex$Win)
    null.pmeans <- pmeans.norm[null.idx, sample.select,drop=F]
    null.pmeans.ecdf <- score.ecdf(null.pmeans, max.winsize=max(sites.js.bfcut.1$size)/winsize+1)
    sites.js.bfcut.1a <- attach.pmeans.score(sites.js.bfcut.1,
                                             index.col,
                                             sites.js.ex,
                                             pmeans.norm,
                                             sample.select)
    sites.js.bfcut.2 <- attach.pmeans.cutoff(sites.js.bfcut.1a,
                                             null.pmeans.ecdf)
    
  } else {
    sites.js.bfcut.2 <- sites.js.bfcut.1
  }
  
  ##bf.ecdf <- e.res.BF.ecdf(e.res, max.winsize=max(sites.js.bfcut.1$size)/winsize+1)

  ##sites.js.bfcut.2 <- attachP.BF.cutoff(sites.js.bfcut.1, bf.ecdf, winsize=winsize)
  ## if(which.score %in% c("score.1","score","score.2")) {
  ##   sites.BF <- subset(sites.js.bfcut.2, subset=sites.js.bfcut.2[,which.score] > bf.cutoff)
  ## } else if(which.score == "score.P") {
  ##   sites.BF <- subset(sites.js.bfcut.2, subset=sites.js.bfcut.2[,which.score] < bf.cutoff)
  ## } else stop("Invalid bf score!")
  sites.BF <- bf.filter.selSites(sites.js.bfcut.2, bf.cutoff, which.score)

  sites.ID <- apply(sites.js.ex[,c("chr","pattern","ID")],1, paste, collapse="-")

  list(sites.js.ex = sites.js.ex[sites.ID%in%sites.BF$siteID,],
       sites.bf = sites.js.bfcut.2)
}
