overlapProp.byID <-
function(dat1,dat2,id1,id2) {
  ##browser()
  stopifnot(nrow(id1)==nrow(id2))
  p1 <- p2 <- rep(NA, nrow(id1))
  na1 <- is.na(id1[,1])
  na2 <- is.na(id2[,1])
  nona <- !(na1 | na2)

  ## subsetting
  ## dat1.idx <- match(paste(id1[,1],id1[,2],sep="."),with(dat1,paste(ID,pattern,sep=".")))
  ## dat2.idx <- match(paste(id2[,1],id2[,2],sep="."),with(dat2,paste(ID,pattern,sep=".")))

  ## dat1.sub <- dat1[dat1.idx,]
  ## dat2.sub <- dat2[dat2.idx,]
  dat.commSub <- overlap.commonSet(dat1, dat2, id1, id2)
  ##id1.nna <- match(id1[nona], dat1$ID)
  ##id2.nna <- match(id2[nona], dat2$ID)

  R1 <- dat.commSub$dat1[nona,c('start','end')]
  R2 <- dat.commSub$dat2[nona,c('start','end')]
  props <- overlap.propotion(R1,R2)

  p1[nona] <- props$p1
  p2[nona] <- props$p2
  data.frame(p1=p1, p2=p2)
}
