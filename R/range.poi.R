range.poi <-
function(s1, s2, which.col=2) {
  if(!missing(s2)) {
    if(nrow(s1) <1 & nrow(s2) <1) return(NA)
    else {
      range.s1 <- if(nrow(s1)==0) NA else range(s1[,which.col])
      range.s2 <- if(nrow(s2)==0) NA else range(s2[,which.col])
      range(range.s1, range.s2, na.rm=T)
    }
  }
  else {
    if(nrow(s1) <1) return(NA)
    else {
      range.s1 <- if(nrow(s1)==0) NA else range(s1[,which.col])
      range.s1
    }
  }
}
