Q3 <-
function(k, m, p) {
  stopifnot(k>2)
  A1 <- 2*b.f(k, m, p)*Fb.f(k-1, m, p) * ((k-1)*Fb.f(k-2,m,p)-m*p*Fb.f(k-3,m-1,p))
  A2 <- .5*b.f(k,m,p)^2*((k-1)*(k-2)*Fb.f(k-3,m,p)-2*(k-2)*m*p*Fb.f(k-4,m-1,p)+m*(m-1)*p^2*Fb.f(k-5,m-2,p))
  A3 <- sum(b.f(2*k-seq(k-1),m,p)*Fb.f(seq(k-1)-1,m,p)^2)
  A4 <- sum(b.f(2*k-seq(k-1),m,p)*b.f(seq(k-1),m,p)*((seq(k-1)-1)*Fb.f(seq(k-1)-2,m,p)-m*p*Fb.f(seq(k-1)-3,m-1,p)))
  Fb.f(k-1,m,p)^3-A1+A2+A3-A4
}
