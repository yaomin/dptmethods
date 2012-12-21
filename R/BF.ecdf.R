BF.ecdf <-
function(max.winsize=200, sample.size=10000, summary.fun=mean, by.chr=F, max.winsize.asym=50) {
  ##browser()
  sample.fun <- function(ae.res,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size, mean=mean(ae.res), sd=sqrt(var(ae.res)/rep.n)))
    }
  }
  chrs <- get.ws.chrs()
  winsize <- get.winsize()
  out.list <- list()
  if(by.chr) {
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      e.res.nrow <- nrow(e.res)
      if(sample.size>e.res.nrow) sample.size <- e.res.nrow
      cat(chr,":",e.res.nrow)
      sfExport("e.res","max.winsize","sample.size","summary.fun","sample.fun")
      out <- sapply(e.res, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
      names(out) <- names(e.res)
      out.list <- c(out.list, list(out))
      cat("\tDONE\n")
    }
    names(out.list) <- chrs
  } else {
    ##browser()
    all.eres <- NULL
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      for(i in seq(e.res)) {
        all.eres[[i]] <- c(all.eres[[i]], e.res[[i]])  
      }
      cat(length(all.eres[[1]]),'\n')
      if(length(all.eres[[1]])>50000) break
      }
    sfExport("all.eres","max.winsize","sample.size","summary.fun","sample.fun")
    out.list <- sapply(all.eres, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
    names(out.list) <- names(e.res)
  }
    
  out.list
}
