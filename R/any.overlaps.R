any.overlaps <-
function(R1, R2, e=500) {
  R1 <- t(apply(R1, 1, `+`, c(-e, e)))
  R2 <- t(apply(R2, 1, `+`, c(-e, e)))
  ret <- NULL
  for ( i in seq(nrow(R1))) {
    ret <- c(ret, any(apply(R2, 1, function(y,x) if.overlap(x,y), x=R1[i,])))
  }
  ret
}
