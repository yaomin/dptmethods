default.pattContrast <-
function(patt) {
  ##browser()
  if(!all(patt==1)) {
    n.patt <- length(patt)
    patt.flip <- rep(1,n.patt)-patt
    if(patt.flip[1]>0) base.contr <- patt
    else base.contr <- patt.flip
    ## if(n.patt>2)
    ##   out.contr <- cbind(base.contr,
    ##                      matrix(0, nrow=n.patt, ncol=n.patt-2))
    ## else out.contr <- as.matrix(base.contr)
    out.contr <- as.matrix(base.contr)
    dimnames(out.contr) <- NULL
   
  }
  else out.contr <- NULL
  out.contr
}
