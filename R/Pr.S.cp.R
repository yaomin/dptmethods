Pr.S.cp <-
function(r,k,n,p, rho=1) {
  exp(-1*n*(r-k*p)*exp(-1*(k*H.f(r/k,p)+rho*h.f(r/k,p)))/sqrt(2*pi*r*k*(k-r)))
}
