compute.qvalues <-
function(pvalues,...) {
  qvalue::qvalue(pvalues, gui=F,...)$qvalues
}
