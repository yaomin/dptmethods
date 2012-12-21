compute.v <-
function(p,p.var,k=1e-36) {
  var.zeros <- min(p.var[p.var !=0])
  p.var[p.var ==0] <- var.zeros
  ##p.zeros <- min(p[p!=0])
  ##p.one <- max(p[p!=1])
  ((p+k)^2*(1-p+k)^2)/((1-p+k)^2+(p+k)^2)*1/p.var
}
