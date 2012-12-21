get.wsoption.para <-
function(which=c("winsize"), pos=-1) {
  which <- match.arg(which)
  if(which == "winsize") get.winsize()
  else stop("Not a valid DiffSeq parameter: ", which)
}
