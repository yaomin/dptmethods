format.data.lme <-
function(dat, group=c(1,1,1,2,2,2,3,3,3), smpid=NULL, log.t=TRUE) {
  y <- as.numeric(unlist(dat))+0
  if(log.t) y <- log2(y)
  win <- rep(row.names(dat),ncol(dat))
  smp <- rep(names(dat), each=nrow(dat))
  grp <- as.factor(rep(group, each=nrow(dat)))
  res <- data.frame(y=y, win=win, smp=smp, grp=grp)
  if(!is.null(smpid)) {
    pid <- as.factor(rep(smpid, each=nrow(dat)))
    res <- cbind(res, smpid=smpid)
  } 
  res
}
