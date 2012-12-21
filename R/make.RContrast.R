make.RContrast <-
function(pattContr) {
  ##browser()
  base.contr <- t(pattContr)[1,]
  zero.idx <- base.contr==0
  one.idx <- base.contr==1
  base.contr[zero.idx]= -1/sum(zero.idx)
  base.contr[one.idx]=1/sum(one.idx)
  out.contr <- make.contrasts(base.contr)
  dimnames(out.contr) <- NULL
  out.contr
}
