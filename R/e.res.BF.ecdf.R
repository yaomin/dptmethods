e.res.BF.ecdf <-
function(e.res, max.winsize=200, sample.size=10000, summary.fun=mean,
                          max.winsize.asym=50, sample.size.asym=1000) {
  ##browser()
  ## sample.fun <- function(ae.res,rep.n) {
  ##   ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
  ## }
  sample.fun <- function(ae.res,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ae.res, min(length(ae.res),sample.size))), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size.asym, mean=mean(ae.res), sd=sqrt(var(ae.res)/rep.n)))
    }
  }
  sapply(e.res, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
}
