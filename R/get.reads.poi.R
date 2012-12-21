get.reads.poi <-
function(reads.poi, start, end) {
  idx.s <- as.numeric(reads.poi) > as.numeric(start)
  idx.e <- as.numeric(reads.poi) < as.numeric(end)
  reads.poi[idx.s & idx.e]
}
