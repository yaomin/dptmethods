#' @export
get.reads <-
function(reads, start, end) {
  idx.s <- as.numeric(reads$poi) >= as.numeric(start)
  idx.e <- as.numeric(reads$poi) <= as.numeric(end)
  reads[idx.s & idx.e,]
}
