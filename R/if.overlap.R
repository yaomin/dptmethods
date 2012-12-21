if.overlap <-
function(r1, r2) {
  ret <- T
  this.order <- order(c(as.numeric(r1), as.numeric(r2)))
  this.order.rev <- order(c(as.numeric(r2), as.numeric(r1)))
  if(sum(this.order-seq(4) != 0) ==0) ret <- F
  else if(sum(this.order.rev-seq(4) != 0) == 0) ret <- F
  ret
}
