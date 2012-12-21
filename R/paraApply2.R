paraApply2 <-
function(x, margin, fun, ..., tmp.dir=".tmp",ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Higher throughput with shared data file
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  tmp <- tempfile(tmpdir=tmp.dir)
  xdim <- dim(x)
  .cut <- cut(seq(xdim[margin]), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  fun.text <- deparse(substitute(fun))
  save(x, margin, xdim, ncore, fun.text, .cut, file=tmp)
 
  .worker <- function(y,load.file,...) {
    ##browser()
    load(load.file)
    fun <- eval(parse(text=fun.text))
    .len <- sum(.cut==levels(.cut)[y])
    if(y>1) {
      .offset <- sum(.cut%in%levels(.cut)[1:(y-1)])
    } else {
      .offset <- 0
    }
    ovec <- apply(x[seq(.len)+.offset,],margin, fun,...)
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    load.file=tmp,
                    simplify=F,
                    ...)
  unlink(tmp)
  unsplit(out.l, .cut)
}
