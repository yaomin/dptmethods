Q2 <-
function(k, m, p) {
  ## (4.4)
  Fb.f(k-1, m, p)^2
  - (k-1)*b.f(k, m, p)*Fb.f(k-2, m, p)
  + m*p*b.f(k, m, p)*Fb.f(k-3, m-1, p)
}
