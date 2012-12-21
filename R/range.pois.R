range.pois <-
function(x, ..., which.col=2) {

  if(nargs()<4) range.poi(x,..., which.col=which.col)
  else {
    range(if(nrow(x) >0 ) range(x[,which.col]) else NA,
          range.pois(..., which.col=which.col), na.rm=T)
  }
}
