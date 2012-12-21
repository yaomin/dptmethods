sel.sigSites <-
function(a.test.report, method="fdr", cutoff=0.05) {
  sigsites <- subset(a.test.report,
                     subset=p.adjust(a.test.report$pvalue, method) < cutoff)
  sig.pos <- subset(sigsites, effect > 0)
  sig.neg <- subset(sigsites, effect < 0)
  list(up=sig.pos, down=sig.neg)
}
