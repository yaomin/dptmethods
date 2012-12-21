add.counts.chr <-
function(chr, x, ...,by=100) {
  arglist <- lapply(list(x,...), function(x) subset(x, subset=(V1 == chr)))
  do.call(add.counts, args=c(arglist, by=by))
}
