#' @export
events.score3 <-
function(res, 
         events.res, 
         rownames,
         sample.label,
	 score.fun="quantile",
	 probs=0.27) {
  ## so called event probability
#  if(length(which.score)>1) stop("Select only one number from 1-5!")
  events.res <- as.data.frame(events.res)
  sel.nonull <- !unlist(lapply(res, is.null)); res <- res[sel.nonull]
  
  p3 <- c(unlist(lapply(res, dptmethods:::get.p)))
  w3.res <- as.matrix(t(laply(res, dptmethods:::get.w3)))
  dimnames(w3.res) <- NULL
  cat("num of input sites:", nrow(w3.res),"\n")
  e.res <- vector('list', ncol(events.res))
  
  e.res <- event.prob3(w3.res, sample.label, events.res, p3, 
		       fun=score.fun,probs=probs,p.fun="median")

  row.names(e.res) <- rownames
  attr(e.res, 'p3') <- p3
  attr(e.res, 'event') <- events.res
  attr(e.res, 'label') <- sample.label
  e.res
}
