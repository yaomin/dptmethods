par.poiComparePerWin <-
function(dat.chr,
                                 dir,
                                 file.pre="",
                                 n.reads,
                                 grp=factor(rep(1:3, each=3)),
                                 n=10000,
                                 n.clusters=2,
                                 sourcefile=NULL)
{
  idx.split <- split(seq(dim(dat.chr)[1]), (seq(dim(dat.chr)[1])-1)%/%n)
  sfInit(parallel = T, cpus = n.clusters)
  sfExport('dat.chr', local=T)
  ##sfSource(sourcefile)
  names(idx.split) <- paste("s", names(idx.split),sep="")
  sfLapply(idx.split,
           function(idx, dat, n.reads, grp, dir, file.pre){
             filenm <- paste(dir,
                             paste(file.pre, paste(range(dat[idx,1]),collapse="-"),"RData", sep="."),
                             sep="/")
             if(file.access(filenm)<0) {
               res <- apply(dat[idx,], 1, poiComparePerWin, n.reads);
               res <- as.data.frame(t(res));
               row.names(res) <- dat[idx,1];
               save(res, file=filenm);
             }
           },
           dat=dat.chr, n.reads=n.reads, grp=grp, dir=dir, file.pre=file.pre)
  sfStop()
  ## res <- as.data.frame(t(as.data.frame(res)))
  ##   row.names(res) <- dat.chr[,1]
  ##   names(res) <- c(sapply(c("est", "std", "t", "P"),
  ##                        function(x, vec) paste(x, vec, sep="."),
  ##                        vec=c("int", "T1vsN","T2vsN","T1vsT2", "T1T2vsN")), "disp", "p.overall")
  ##   res
}
