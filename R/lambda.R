lambda <-
function(r,k,n,p) {
  res <- (n-k+1)*(r/k-p)*dbinom(r,k,p)
}
