compute.minB.p <-
function(x, paras) {
  a <- paras$ab[,1]
  b <- paras$ab[,2]
  e <- paras$e
  q <- (2^x) * e/(1 - e + (2^x)*e)
  pb <- pbeta(q, a, b)
  ##print(pb)
  prod(1-pb)
}
