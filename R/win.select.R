#' @export
win.select <-
function(e.res, e.cut, exclude.zeroPattern=T, zeroPattern.cutoff=0){
  ##browser()
  method = attr(e.cut, 'method')
  stopifnot(length(e.res)==length(e.cut))
  ## select the windows by thresholding on event probabilities or BF
  p3 <- attr(e.res, 'p3')
  para <- attr(e.cut, "para")
  events.res <- attr(e.res, 'event')
  n.win <- nrow(e.res)
  n.event <- ncol(e.res)
  wins <- row.names(e.res)

  wins.sel <- vector('list', n.event)
  for (i in seq(n.event)) {
    sel.i <- e.res[,i]> e.cut[i]
    scores.i <- e.res[sel.i,i]
    wins.i <- as.numeric(wins[sel.i])
    p.i <- vector('numeric', length(scores.i))
    if(method=="beta" & length(scores.i)>0) p.i <- sapply(scores.i,compute.minB.p, para[[i]])
    else {
      ecdf.i <- ecdf(e.res[,i])
      p.i <- 1-ecdf.i(scores.i)
    }
    wins.sel[[i]] <- data.frame(scores=scores.i, wins=wins.i, probs=p.i,row.names=wins.i)    
  }
  names(wins.sel) <- names(e.res)
  attr(wins.sel, 'event') <- events.res
  attr(wins.sel, 'cutoff') <- e.cut
  ## removing zero patterns from other if asked
  if(exclude.zeroPattern) {
    idx.0 <- which.zeroPattern(names(e.res))
    idx.oth <- setdiff(seq(n.event), idx.0)
    wins.zero <- wins[e.res[,idx.0]>zeroPattern.cutoff]
    for (j in idx.oth) {
      wins.sel[[j]] <- subset(wins.sel[[j]],subset=!(wins.sel[[j]]$wins%in%wins.zero))
    }
  }
  wins.sel
}
