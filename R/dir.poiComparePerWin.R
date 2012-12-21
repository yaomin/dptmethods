dir.poiComparePerWin <-
function(dir,
                                 n.reads,
                                 grp=factor(rep(1:3, each=3)),
                                 n=100000,
                                 n.clusters=2) {
  files <- grep("JoinedSamples.chr", list.files("../data/Colon"), v=T)
  for (eachFile in files) {
    datnm <- load(paste(dir, eachFile, sep="/"))
    dat <- consolidate.lows(get(datnm), cols=2:10)
    if(dim(dat)[1]>0) {
      res <- par.poiComparePerWin(dat,
                                  dir=dir,
                                  file.pre=datnm,
                                  n.reads=n.reads,
                                  n=n,
                                  n.clusters=n.clusters)
      ## save(res, file=paste(dir, paste("poi",datnm, "RData", sep="."), sep="/"))
      ## cat(paste(datnm,"is done","\n"))
    }
  }
}
