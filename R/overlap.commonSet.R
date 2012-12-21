overlap.commonSet <-
function(dat1, dat2, id1, id2) {
  stopifnot(nrow(id1)==nrow(id2))
  ## subsetting
  dat1.idx <- match(paste(id1[,1],id1[,2],sep="."),with(dat1,paste(ID,pattern,sep=".")))
  dat2.idx <- match(paste(id2[,1],id2[,2],sep="."),with(dat2,paste(ID,pattern,sep=".")))

  dat1.sub <- dat1[dat1.idx,]
  dat2.sub <- dat2[dat2.idx,]

  list(dat1=dat1.sub, dat2=dat2.sub)
}
