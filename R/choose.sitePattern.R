choose.sitePattern <-
function(sites, e.cutoff=-1, site.cutoff=0.7) {
  ##browser()
  chrs <- as.character(unique(sites$chr))
    ##print(chrs[i])
  .out <- sfSapply(chrs,
                   choose.sitePattern.byChr,
                   sites=sites,
                   e.cutoff=e.cutoff,
                   site.cutoff=site.cutoff,
                   simplify=F)
  ##ldply(.out, function(x) x)
  .out
}
