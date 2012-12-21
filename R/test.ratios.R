test.ratios <-
function(dat, x.n, y.n, w.n=2500000){
 ## Assume we have data from join.samples
 ##require(rateratio.test)

 res <- apply(dat,
              1,
              function(adat, x.n,y.n, w.n)
              as.numeric(unlist(rateratio.test(adat, c(x.n, y.n)/w.n))[c(1,2,3,4,6,7)]),
              x.n=x.n,
              y.n=y.n,
              w.n=w.n)
 t(res)
}
