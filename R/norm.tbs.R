norm.tbs <-
function(tbs, bg.correct=F) {
  ##browser()
  tb.dn <- denoiseSignal(tbs)
  bg <- bg.N(tbs)
  ## tb.norm <- lapply(tb.dn, function(x) {x$freq <- x$freq + bg[x$count+1,'freq'];
  ##                                       x=rbind(c(0, bg[bg$count==0, 'freq']),x);
  ##                                     x})
  if(bg.correct) {
    tb.norm <- tb.dn
  } else {
    tb.norm <- lapply(tb.dn, function(x) {x$freq <- x$freq + bg[x$count+1,'freq'];x})
  }
   
  tb.norm
}
