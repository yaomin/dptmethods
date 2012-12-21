assign.groups <-
function(pmeans, cuts, patt=c(0,0,0,1,1,1,1,1,1)){
  res <- pmeans
  for (i in seq(ncol(pmeans))) {
    res[,i] <- if(patt[i] ==0) pmeans[,i]<= cuts$lo[i] else pmeans[,i] > cuts$hi[i]
  }
  apply(res,1, all)
}
