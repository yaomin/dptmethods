compute.beta.para.3 <-
function(res, rev=F) {
  if(rev) {
    ret <- lapply(res, function(x) if(!is.null(x)) est.beta.para(1-x$w.win[,3]))
  }
  else {
    ret <- lapply(res, function(x) if(!is.null(x)) est.beta.para(x$w.win[,3]))
  }
  ret <- as.data.frame(matrix(unlist(ret), byrow=T, ncol=2))
  names(ret) <- c("a", "b")
  ret
}
