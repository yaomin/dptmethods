poiComparePerWin <-
function(y, n.reads,grp=factor(rep(1:3, each=3))) {
 
  grp1 <- grp2 <- grp
  contrasts(grp1) <- cbind(c(0,1,0),
                           c(0,0,1))

  contrasts(grp2) <- cbind(c(0,1,1),
                           c(0,-1,1))
  dat1 <- data.frame(y=unlist(y[-1]), grp=grp1)
  dat2 <- data.frame(y=unlist(y[-1]), grp=grp2)

  ctl.idx <- which(grp==1)
  tst.ct1 <- glm(y~grp, data=dat1, family=quasipoisson,
                 offset=log(n.reads)-log(mean(n.reads[ctl.idx])))
  tst.ct2 <- glm(y~grp, data=dat2, family=quasipoisson,
                 offset=log(n.reads)-log(mean(n.reads[ctl.idx])))
  res1 <- as.vector(summary(tst.ct1)$coefficients)
  res2 <- as.vector(summary(tst.ct2)$coefficients[-1,])
  disp <- as.vector(summary(tst.ct1)$dispersion)
  overallP <- drop1(tst.ct1, test="F")$'Pr(F)'[2]
  res <- c(res1[1:3],res2[1:2], # estimate
           res1[4:6], res2[3:4], # std
           res1[7:9], res2[5:6], # t value
           res1[10:12], res2[7:8], # Pr
           disp, # dispersion
           overallP) # F test on grp
  res
}
