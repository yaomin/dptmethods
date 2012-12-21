dpt.syshold <-
function(tm, verbose=F) {
  t1 <- proc.time()
  t.lag <- 0
  while(t.lag<tm) {
    system.time(for(i in 1:10000) x <- mean(rt(1000, df=4)))
    t2 <- proc.time()
    t.lag <- (t2-t1)[3]
  }
  if(verbose) cat("system waited for",t.lag,"sec\n")
}
