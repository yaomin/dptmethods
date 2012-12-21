#' @export
get.mappedCountTab <-
function(inPath=NULL, pos=-1) {
  if(is.null(inPath)) {
    inPath <- get.wsoption.path("joinSample", pos=pos)
  }
  cntab.data.file <- file.path(inPath, "countdistr.RData")
  load(cntab.data.file)
  ##assign("counts.tab", count.tbs, pos=parent.frame(n=1))
  lapply(count.tbs, function(x) {x$count <- as.numeric(as.character(x$count));x})
  ##count.tbs
}
