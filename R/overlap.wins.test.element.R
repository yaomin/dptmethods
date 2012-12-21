overlap.wins.test.element <-
function(x1, x2, e, self=F) {
  poi.cols <- c('start','end')
  output <- any.overlap.sites(x1[,poi.cols],x2[,poi.cols],e=e,self=self)
  output
}
