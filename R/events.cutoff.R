events.cutoff <-
function(res, e.res, cutoff, method = c("beta", "BF", "best"),find.range= c(-20,20)) {
  ## Beta prior based distribution
  ## events =
  ##   0 0
  ##   1 0
  ##   1 1
  ## events are defined by columns
  ##browser()

  method <- match.arg(method)
  events <- as.data.frame(attr(e.res, 'event'))
  sample.label <- attr(e.res, 'label')
  if(method=="beta") {
    res <- res[!unlist(lapply(res, is.null))]
    para.names <- paste("P",apply(events, 2, paste, collapse=""), sep="")
    events <- as.data.frame(pattern.design(events, sample.label))
    para <-  apply(events,
                   2,
                   function(x) compute.minB.para(res, x)
                   )
    names(para) <- para.names
    
    e.cut <- lapply(para,
                    function(x) find.cutoff(cutoff,
                                            x,
                                            find.range)
                    )
    
    ret <- unlist(e.cut)
    attr(ret, 'para') <- para
    names(ret) <- NULL
  }
  else if(method=="BF") {
    ret <- rep(cutoff, ncol(events))
    attr(ret, 'para') <- NA
  }
  else if(method=="best") {
    .cutoff <- findBest.escoreCutoff(e.res, seq(-2,2,0.01), cutoff)
    ret <- rep(.cutoff, ncol(events))
    attr(ret, 'para') <- NA
  }
  attr(ret, 'method') <- method
  ret
}
