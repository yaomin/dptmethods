join.counts.chr <-
function(chr, x, ...,by=100) {
  
  ## notice: the data structures of x,... are differet from those in join.samples.chr
  arglist <- lapply(list(x,...), function(x) x[[chr]])
  do.call(join.counts, args=c(arglist,by=by))
}
