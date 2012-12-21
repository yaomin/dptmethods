BF.cutoff <-
function(max.winsize=200, sample.size=10000, summary.fun=mean, cutoff=0.001, by.chr=F) {
  ##browser()
  sample.fun <- function(ae.res,rep.n) {
    quantile(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun), probs=1-cutoff)
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
      sfExport("e.res","max.winsize","sample.size","summary.fun","cutoff","sample.fun")
      ##sfLibrary(snowfall)
      out <- lapply(e.res, function(x) {##sfExport("x","max.winsize","sample.size","summary.fun","cutoff","sample.fun");
        sfSapply(seq(max.winsize), function(y) sample.fun(x,y))})
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
    sfExport("all.eres","max.winsize","sample.size","summary.fun","cutoff","sample.fun")
    out.list <- lapply(all.eres, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)))
    ##browser()
    names(out.list) <- names(e.res)
  }
  out.1 <- sapply(out.list,
                  function(x) ldply(x, function(y) data.frame(y)),
                  simplify=F)
  out.df <- ldply(out.1,
                  function(z) {names(z) <- c('pattern','x');
                               z}
                  )
  if(by.chr) names(out.df) <- c("chr","pattern","cutoff")
  else {out.df <- out.df[,c(1,3)]
        names(out.df) <- c("pattern", "cutoff")
      }
  out.df
}
