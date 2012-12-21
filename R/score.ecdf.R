score.ecdf <-
function(score, max.winsize=200, sample.size=10000, summary.fun=mean,
                          max.winsize.asym=50, sample.size.asym=1000) {
  ##browser()
  ## sample.fun <- function(ae.res,rep.n) {
  ##   ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
  ## }
  sample.fun <- function(ascore,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ascore, min(length(ascore),sample.size))), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size.asym, mean=mean(ascore), sd=sqrt(var(ascore)/rep.n)))
    }
  }
  sapply(score, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
}
