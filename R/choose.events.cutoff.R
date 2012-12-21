#'  @export

choose.events.cutoff <-
function(res, e.res, cutoff, method=c("BF","beta","best")) {
  method <- match.arg(method)
  ## if(method == "BF") out <- events.cutoff(res, e.res, cutoff, 'BF')
  ## else if(method == "beta") out <- events.cutoff(res, e.res, cutoff, 'beta')
  out <- events.cutoff(res, e.res, cutoff, method)
  out
}
