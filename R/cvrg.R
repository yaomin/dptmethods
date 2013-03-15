## Filename: cvrg-fun.r
## By: Yaomin Xu (yaomin@gmail.com)
## Date: Thu Jan 10 16:52:18 2013
## Notes:   Functions called in the coverage workflow
## =============================================================================

bamHeaders <- function(files) {
  goodH <- list()
  goodF <- NULL
  for (i in seq(along=files)) {
    aHeader <- tryCatch(scanBamHeader(files[i]),
                        error=function(e) e)
    if(!inherits(aHeader, "error")) {
      goodH <- c(goodH, aHeader)
      goodF <- c(goodF, files[i])
    }
  }
  list(files=goodF, headers=goodH)
}

para.countBam <- function(files, ncore=3,...) {
  cl <- makeCluster(getOption("cl.cores",ncore), type="FORK")
  cnts <- parLapplyLB(cl,
                      files,
                      Rsamtools::countBam,
                      ...)
  stopCluster(cl)
  names(cnts) <- as.character(sapply(cnts, function(x) x[['file']]))
  if(length(files)!=length(names(cnts))) warning("Not all input files read!")
  cnts
}

scanBamFiles <- function(files, ncore=3,...) {
  bfs.headers <- bamHeaders(files)
  good.bfs <- bfs.headers$files
  chrlens <- bfs.headers$headers[[1]]$target
  bfs.counts <- para.countBam(good.bfs, ncore=ncore)
  out <- list(headers=bfs.headers,
              files=good.bfs,
              chrlens=chrlens,
              recordCount=bfs.counts)
  out
}

bam2coverage <- function(file,
                         chrlens=NULL,
                         ext=0L,
                         coords=c("leftmost", "fiveprime"),
                         format="BAM",
                         flag=NULL,
                         shift=0L)
{
  coords <- match.arg(coords)

  if(is.null(flag)) {
    flag <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
  }

  scanParam <- ScanBamParam(flag=flag)
  
  aligns <- readGappedAlignments(file, use.names=F, param=scanParam)

  bpcv <- dpt.bpcoverage(aligns,
                         ext=ext,
                         w=chrlens,
                         shift=shift,
                         uniq=F,
                         coords=coords)
  list(data=bpcv,
       metadata=list(length=length(aligns), seqinfo=seqinfo(aligns)))
}

para.bam2coverage <- function(files,
                              ncore=3,
                              use.names=T,
                              chrlens=NULL,
                              ext=0L,
                              coords=c("leftmost", "fiveprime"),
                              format="BAM",
                              flag=NULL,
                              shift=0L)
{

  cl <- makeCluster(getOption("cl.cores",ncore), type="FORK")
  cvl <- parLapplyLB(cl,
                     files,
                     bam2coverage,
                     chrlens=chrlens,
                     ext=ext,
                     coords=coords,
                     format=format,
                     flag=flag,
                     shift=shift)
  stopCluster(cl)
  if(use.names) names(cvl) <- basename(files)
  cvl
}

dpt.bpcoverage <- function(x, shift=0L, w=NULL, ext=0L, uniq=F, coords=c("leftmost","fiveprime"),...) {
  ## x is a GappedAlignments
  ## w is chromosome width with names, it also filters chrs
  
  coords <- match.arg(coords)
  rgs <- ranges(x)
  strnd <- as.vector(strand(x))
  if(is.null(w)) w <- seqlengths(x)
  
  if(coords=="leftmost") {
    rstart <- start(rgs)-(strnd=="-")*ext
    rend <- end(rgs)+(strnd=="+")*ext
  } else {
    rstart <- start(rgs)-(strand(x)=="-")*(width(rgs)+ext-1L)
    rend <- start(x)+(strand(x)=="+")*(width(rgs)+ext-1L)
  }
  start(rgs) <- rstart
  end(rgs) <- rend
  rls <- split(rgs, seqnames(x))

  if(is.null(w)) {
    chrs <- levels(seqnames(x))
  } else {
    chrs <- names(w)
  }                 
  
  cvg <- RleList(lapply(chrs,
                        function(chr, ...) {
                          chrIR <- rls[[chr]]
                          if (identical(shift, 0L)) {
                            chr_shift <- 0L
                          } else {
                            chr_shift <- shift[chr]
                          }
                          if (is.null(w)) {
                            chr_width <- NULL
                          } else  {
                            chr_width <- w[chr]
                          }
                          if (uniq) {
                            IRanges::coverage(unique(chrIR),
                                              shift=chr_shift, width=chr_width, ...)
                          } else {
                            IRanges::coverage(chrIR,
                                              shift=chr_shift, width=chr_width, ...)
                          } 
                        },
                        ...),
                 compress = FALSE)
  names(cvg) <- chrs
  cvg
}

get.data.bfscvrg <- function(cvrg) {
  counts <- lapply(cvrg, function(x) x$data)
  counts
}


get.metadata.bfscvrg <- function(cvrg, what=c("length","seqinfo")) {
  what <- match.arg(what)
  lens <- lapply(cvrg, function(x) x$metadata[[what]])
  lens
}

  
