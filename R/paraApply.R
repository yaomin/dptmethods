paraApply <-
function(x, margin, fun, ...,mode=real64(),mode.size=8, tmp.dir=".tmp",ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Memory efficient with mmap
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  if(missing(mode.size)) mode.size <- sizeof(mode)
  tmp <- tempfile(tmpdir=tmp.dir)
  xdim <- dim(x)
  .cut <- cut(seq(xdim[margin]), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  
  if(margin==1) {
    writeBin(as.vector(t(x)), tmp)
    frac.len <- ncol(x)
  } else {
    writeBin(as.vector(x), tmp)
    frac.len <- nrow(x)
  }
  
  fun.text <- deparse(substitute(fun))
  ##in.mmap <- mmap(tmp, mode=mode, )

  .worker <- function(y,
                      margin,
                      .file.mmap,
                      frac.len,
                      mode,
                      mode.size,
                      xdim,
                      ncore,
                      fun.text,
                      .cut,
                      ...) {
    ##browser()
    fun <- eval(parse(text=fun.text))
    .len <- sum(.cut==levels(.cut)[y])
    if(y>1) {
      .offset <- sum(.cut%in%levels(.cut)[1:(y-1)])
    } else {
      .offset <- 0
    }
    .data.mmap <- mmap(.file.mmap, mode=mode, prot=mmapFlags("PROT_READ"))
    ##cat(class(length(.data.mmap)))
    ##cat("which=",which)
    ovec <- tapply(.data.mmap[seq(.len*frac.len)+.offset*frac.len],
                   list(rep(seq(.len), each=frac.len)), fun,...)
    dimnames(ovec) <- NULL
    munmap(.data.mmap)
    ##cat(str(ovec))
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    margin=margin,
                    .file.mmap=tmp,
                    frac.len=frac.len,
                    mode=mode,
                    mode.size=mode.size,
                    xdim=xdim,
                    ncore=ncore,
                    fun.text=fun.text,
                    .cut=.cut,
                    simplify=F,
                    ...)
  ##munmap(in.mmap)
  unlink(tmp)
  unsplit(out.l, .cut)
}
