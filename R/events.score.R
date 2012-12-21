#' @export
events.score <-
function(res, 
                         events.res, 
                         rownames, 
                         sample.label, 
                         which.score=2, 
                         npart=10, 
                         para.mode=2) {
  ## so called event probability

  if(length(which.score)>1) stop("Select only one number from 1-5!")
  events.res <- as.data.frame(events.res)
  sel.nonull <- !unlist(lapply(res, is.null)); res <- res[sel.nonull]
  
  p3 <- c(unlist(lapply(res, get.p)))
  w3.res <- as.matrix(t(laply(res, get.w3)))
  dimnames(w3.res) <- NULL
  cat("num of input sites:", nrow(w3.res),"\n")
  e.res <- vector('list', ncol(events.res))
 
  patts <- pattern.design(events.res, sample.label)
  for (i in seq(ncol(events.res))) {
    cat("pattern:",paste(patts[,i],sep="",collapse=""))
    .this <- event.prob(w3.res, 
                        sample.label,
                        event=patts[,i],
                        p=p3, 
                        which.score=which.score, 
                        npart=npart,
                        para.mode=para.mode)
    e.res[[i]] <- as.vector(.this)
    cat("\t","done\n")
  }
  e.res <- as.data.frame(e.res)
  names(e.res) <- paste("P", apply(events.res,2,paste, collapse=""), sep="")
  row.names(e.res) <- rownames
  attr(e.res, 'p3') <- p3
  attr(e.res, 'event') <- events.res
  attr(e.res, 'label') <- sample.label
  e.res
}
