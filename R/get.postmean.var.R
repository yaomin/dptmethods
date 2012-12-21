get.postmean.var <-
function(x) {
  x$signal.var/(x$signal.win)^2
}
