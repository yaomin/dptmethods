## Filename: binning-fun.r
## By: Yaomin Xu (yaomin@gmail.com)
## Date: Thu Jan 10 17:13:00 2013
## Notes:   Functions called in binning process
## =============================================================================
dpt.binCoverage <- function(cvg, bsize=50) {
  bcvg <- RleList(lapply(cvg, binMeanVector, bsize=bsize), compress=FALSE)
  names(bcvg) <- names(cvg)
  bcvg
}

binMeanVector <- function (vec, bsize) 
{
  m <- length(vec)%/%bsize
  r <- length(vec)%%bsize
  ##   starts <- bsize * 0:m + 1
  ##   ends   <- c(bsize * 1:m,length(vec))
  mean1 <- colMeans(matrix(as.vector(vec)[seq(m*bsize)], nrow=bsize))
  mean2 <- mean(as.vector(vec)[(m*bsize):length(vec)])
  ##RangedData(IRanges(start=starts, end=ends), space=chr, coverage=c(mean1,mean2))
  c(mean1, mean2)
}

para.binCoverage <- function(cl, cvrg, ncore, bsize) {
##  cl <- makeCluster(getOption("cl.cores",ncore), type="FORK")
  bcvrg <- parLapplyLB(cl,
                       cvrg,
                       dpt.binCoverage,
                       bsize=bsize)
##  stopCluster(cl)
  list(data=bcvrg, metadata=list(bsize=bsize))
}

get.bcvrg <- function(x, what) {
  x[[what]]
}
  
get.chr.bcvrg <- function(bcvrg, chr, returnType=c("DataFrame","list")) {
  returnType <- match.arg(returnType)
  out <- lapply(get.bcvrg(bcvrg, "data"), function(x) x[[chr]])
  if(returnType=="DataFrame") {
    as(out, "DataFrame")
  } else out
}

round.bcvrg <- function(x) {
  .data <- lapply(get.bcvrg(x, "data"),
                  function(y) lapply(y, function(z) round(z,0)))
  list(data=.data, metadata=get.bcvrg(x, "metadata"))
}

fpkm.bcvrg <- function(x, lens,k=1e3, m=1e6) {

  .dat <- get.bcvrg(x, "data")
  .out <- lapply(names(lens), function(.smp) lapply(.dat[[.smp]], function(y) y*k*m/lens[[.smp]]))
  ## .out <- list()
  ## for(.smp in names(lens)) {
  ##   .out <- list(.out, lapply(.dat[[.smp]], function(y) y*k*m/lens[[.smp]]))
  ## }
  names(.out) <- names(lens)
  list(data=.out, metadata=get.bcvrg(x, "metadata"))
}

bcvrg.DataFrameList <- function(bcvrg, chrs=NULL) {
  if(is.null(chrs)) chrs <- names(bcvrg$data[[1]])
  data <- DataFrameList(lapply(chrs, get.chr.bcvrg, bcvrg=bcvrg))
  names(data) <- chrs
  list(data=data, metadata=bcvrg$metadata)
}


len2IRanges <- function(len, bsize) {
  m <- len%/%bsize
  starts <- bsize * 0:m + 1
  ends   <- c(bsize * 1:m,len)
  IRanges(start=starts, end=ends)
}

get.bcvrg.DFList <- function(bcvrg.dfl, what=c("data","metadata")) {
  what <- match.arg(what)
  bcvrg.dfl[[what]]
}

bcvrg.DFList2RD <- function(bcvrg.dfl, chrlens) {
  ##browser()
  .data <- get.bcvrg.DFList(bcvrg.dfl, "data")
  bsize <- get.bcvrg.DFList(bcvrg.dfl,"metadata")[['bsize']]
  comchrs <- intersect(names(.data), names(chrlens))
  .rlist <- RangesList(lapply(chrlens[comchrs], len2IRanges, bsize))
  .out <- RangedData(.rlist, get.bcvrg.DFList(bcvrg.dfl, "data")[comchrs])
  .out
}

colApply.bcvrg.dfl <- function(x, fun, ...) {
  .values <- DataFrameList(lapply(values(x),
                                  function(y) apply(as.data.frame(y),
                                                    2,
                                                    fun,...)))
  .rlist <- ranges(x)
  RangedData(.rlist, .values)
}

cbind.bcvrg.dfl <- function(...) {
  .list <- lapply(list(...), get.bcvrg.DFList, what="data")
  .dat <- .list[[1]]
  for(alst in .list[-1]) {
    .dat <- cbind(.dat, alst)
  }
  .metadata <- get.bcvrg.DFList(..1, what="metadata")
  list(data=.dat, metadata=.metadata)
}

vnames.bcvrg.dfl <- function(x, vnames) {
  x$data <- DataFrameList(lapply(x$data,
                                 function(y) {
                                   names(y) <- vnames
                                   y
                                 }))
  x
}

process.bcvrg <- function(bcvrg, bsize, ncore=3) {
  ## bcvrg could be bcvrg or bcvrg.fpkm
  
#  bcvrg <- para.binCoverage(cvrg, ncore, bin.size)
  bcvrg.round <- round.bcvrg(bcvrg)
#  bcvrg.fpkm <- fpkm.bcvrg(bcvrg, get.metadata.bfscvrg(bfs.cvrg),1000, 1e6)
  #bcvrg.fpkm.round <- round.bcvrg(fpkm)

  bcvrg.dfl <- bcvrg.DataFrameList(bcvrg)
  bcvrg.dfl.round <- bcvrg.DataFrameList(bcvrg.round)
  #bcvrg.dfl.fpkm <- bcvrg.DataFrameList(fpkm)
  #bcvrg.dfl.fpkm.round <- bcvrg.DataFrameList(bcvrg.fpkm.round)
  
  bcvrg.full <- cbind.bcvrg.dfl(bcvrg.dfl.round,
                                bcvrg.dfl
                                #bcvrg.dfl.fpkm.round,
                                #bcvrg.dfl.fpkm
                                )
  ncols <- ncol(bcvrg.full[[1]][[1]])
  bcvrg.full <- vnames.bcvrg.dfl(bcvrg.full,
                                 paste(##rep(c("s","l","fs","fl"),
                                       rep(c("s","l"),
                                           each=ncols/2),
                                       rep(1:(ncols/2), 2),
                                       sep=""))
  bcvrg.full
}

io.joinsample <- function(bcvrg.full.rd,
                          dir=".",
                          select=NULL,
                          chrs=NULL,
                          offset=0,
                          prefix="joinsample") {
  ##browser()
  if(is.null(select)) select=grep("^s",
                        colnames(bcvrg.full.rd),
                        v=T)
  if(is.null(chrs)) chrs <- names(bcvrg.full.rd)
  renames <- c("poi", select[-1])
  chrs <- names(bcvrg.full.rd)
  cat("Output joinsample dataset:\n")
  for(chr in chrs) {
    cat(chr, "\t")
    .this.rd <- bcvrg.full.rd[chr][,select]
    poi <- if(offset==0) start(.this.rd)-1L else start(.this.rd)
    .this <- data.frame(poi, as.data.frame(.this.rd)[,select])
    .vname <- paste("s", chr, sep=".")
    .fname <- paste(prefix, chr, "RData", sep=".")
    assign(.vname, .this)
    save(list=.vname, file=file.path(dir,.fname))
    cat("done!\n")
  }
}
           
         
                   
