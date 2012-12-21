QL <-
function(k,m,L,p) {
  ## N=m*L
  exp(log(Q2(k,m,p))+(L-2)*(log(Q3(k,m,p))-log(Q2(k,m,p))))
}
