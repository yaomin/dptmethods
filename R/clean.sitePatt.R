#' @export
clean.sitePatt <-
function(patt, sites, e.TF, cutoff=0.5, debug=F) {
  ### sites is IRanges
  ### e.TF is RangedData from ranged.e.TF
  if(debug) browser()
  sites <- sites[[patt]]
  e.sub <- e.TF[ranges(e.TF)[[1]] %over% sites,]
  .fac <- rep(NA, length(e.sub[[1]]))
  e.ranges <- ranges(e.sub)[[1]]
  .fac <- findOverlaps(e.ranges, sites, select="first")
  .summary <-  by(as.data.frame(values(e.sub))[,-1],
                  list(.fac),
                  function(x, cutoff) patt.summary(x, cutoff),
                  cutoff=cutoff,
                  simplify=F)
  .summary.df <- pattSummary.2dataframe(.summary)
  idx.included <- with(.summary.df,
                       best.pattern==patt & best.valid)
  sites.rd <- RangedData(sites, score=.summary.df$diff.score*log(.summary.df$N))
  sites.rd[idx.included,]
}
