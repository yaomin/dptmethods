normalize.signal <-
function(x, tb, tb.norm, out.int=F) {
  ##browser()
  ls0 <- approxfun(tb$count, tb$freq)
  lsnorm <- approxfun(tb.norm$count,tb.norm$freq)

  out <- x*lsnorm(x)/ls0(x)
  if(out.int) out <- as.integer(round(out, 0))
  out
}
