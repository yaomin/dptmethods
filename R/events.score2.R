#' @export
events.score2 <-
function(res, 
         events.res, 
         rownames,
         sample.label,
	 score.fun="mean") {
  ## so called event probability
#  if(length(which.score)>1) stop("Select only one number from 1-5!")
  events.res <- as.data.frame(events.res)
  sel.nonull <- !unlist(lapply(res, is.null)); res <- res[sel.nonull]
  
  p3 <- c(unlist(lapply(res, get.p)))
  w3.res <- as.matrix(t(laply(res, get.w3)))
  dimnames(w3.res) <- NULL
  cat("num of input sites:", nrow(w3.res),"\n")
  e.res <- vector('list', ncol(events.res))
  
  e.res <- event.prob2(w3.res, sample.label, events.res, p3, score.fun)

  row.names(e.res) <- rownames
  attr(e.res, 'p3') <- p3
  attr(e.res, 'event') <- events.res
  attr(e.res, 'label') <- sample.label
  e.res
}
