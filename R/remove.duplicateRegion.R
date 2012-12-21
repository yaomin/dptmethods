remove.duplicateRegion <-
function(x) {
  whichmin <- which.min(x$rank)
  cbind(x[whichmin,], nrow(x))
}
