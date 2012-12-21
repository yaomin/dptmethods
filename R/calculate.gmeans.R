
#'  @export

calculate.gmeans <-
function(tst.regions) {
  
  mx01 <-  ldply(strsplit(as.character(tst.regions$contrast), split=character()), function(x) x)
  n0 <- apply(mx01, 1, function(x) sum(x=='0'))
  n1 <- apply(mx01, 1, function(x) sum(x=='1'))
  n <- n0+n1
  y1 <- tst.regions$value.0+n0/n*tst.regions$value.1
  y0 <- y1 - tst.regions$value.1
  mean0 <- y0
  mean1 <- y1
  ##mean0[tst.regions$value.1>0] <- y0[tst.regions$value.1>0]
  mean0[tst.regions$value.1<=0] <- y1[tst.regions$value.1<=0]
  ##mean1[tst.regions$value.1>0] <- y1[tst.regions$value.1>0]
  mean1[tst.regions$value.1<=0] <- y0[tst.regions$value.1<=0]

  data.frame(mean.0=mean0,
             mean.1=mean1)
}
