#' @export
event.patternsToScan <-
function(allpatt, exclude=NULL,noZeroPatt=T) {
  
  out <- allpatt
  if(!is.null(exclude)) out <- out[,-1*exclude,drop=F]
  if(noZeroPatt & length(which(colSums(out) == 0))>0) out <- out[,-which(colSums(out) == 0),drop=F]
  out
}
