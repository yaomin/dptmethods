clean.sitePatt <-
function(patt, sites, e.TF, cutoff=0.5) {
  ### sites is IRanges
  ### e.TF is RangedData from ranged.e.TF
  ##browser()
  sites <- sites[[patt]]
  e.sub <- e.TF[ranges(e.TF)[[1]] %in% sites,]
  .fac <- rep(NA, length(e.sub[[1]]))
  e.ranges <- ranges(e.sub)[[1]]
  .fac <- match(e.ranges, sites)
  ## useful but not very efficient in this case.
  ##e.TF.rd <- RangedData(ranges=e.ranges, values(e.sub), space=.fac)
  ##.RDPara <- RDApplyParams(e.TF.rd, patt.summary)
  ##.summ <- rdapply(.RDPara)
  .summary <-  by(as.data.frame(values(e.sub))[,-1],
                  list(.fac),
                  function(x, cutoff) patt.summary(x, cutoff),
                  cutoff=cutoff,
                  simplify=F)
  .summary.df <- pattSummary.2dataframe(.summary)
  idx.included <- with(.summary.df,
                       best.pattern==patt & best.valid)
  sites.rd <- RangedData(sites, score=.summary.df$diff.score*log(.summary.df$N))
  ## weighted by width to avoid small site noise
  sites.rd[idx.included,]
}
