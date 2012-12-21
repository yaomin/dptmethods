#' @export
event.allPossiblePatterns <-
function(smplabel) {
  alabel <- smplabel
  grp.n <- length(unique(alabel))
  as.data.frame(t(expand.grid(as.data.frame(matrix(rep(c(0,1),grp.n), ncol=grp.n)))))
}
