f.f <-
function(r,k,p) {
  p/r*dbinom(r-1,k-1,p)*(r*(1-p)*dbinom(r-1,k-1,p)+(r-k*p)*pbinom(r-2,k-1,p))
}
