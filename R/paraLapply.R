paraLapply <-
function(x, fun, ..., tmp.dir=".tmp", ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Higher throughput with shared data file
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  tmp <- tempfile(tmpdir=tmp.dir)
  xn <- length(x)
  .cut <- cut(seq(xn), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  fun.text <- deparse(substitute(fun))
  save(x, xn, ncore, fun.text, .cut, file=tmp)
 
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
    ovec <- lapply(x[seq(.len)+.offset,], fun,...)
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    load.file=tmp,
                    ...)
  unlink(tmp)
  out.l
}
