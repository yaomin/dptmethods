pattern.design <-
function(pattern, sample.label) {
  sl <- factor(sample.label)
  patt <- as.data.frame(pattern)
  stopifnot(nlevels(sl) == nrow(patt))
  levels(sl) <- seq(nlevels(sl))
  ret <- apply(patt, 2, function(x) x[sl])
  if(is.null(dim(ret))) ret <- matrix(ret, nrow=1)
  dimnames(ret) <- NULL
  ret
}
