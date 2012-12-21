choose.sitePattern.byChr <-
function(chrom, sites, e.cutoff=-1, site.cutoff=0.66) {
  ##browser()
  get.ws(ws="pattRecog",
         what="e.res",
         chr=chrom)
  .sites <- subset(sites, subset=chr==chrom)
  ## .out <- apply(.sites, 1, function(x) {cat(x['.id'],"\t");
  ##                                       cat(x['start'], x['end'], '\n');
  ##                                       choose.site.pattern(e.res,x['start'], x['end'],
  ##                                                           e.cutoff,
  ##                                                           site.cutoff)})
  .out <- apply(.sites, 1, function(x) {choose.site.pattern(e.res, x['start'], x['end'],
                                                            e.cutoff,
                                                            site.cutoff)})

  data.frame(ldply(.out, function(x) x$counts),
             N=ldply(.out, function(x) x$N)[,2],
             best.cutoff=ldply(.out, function(x) x$best.cutoff)[,2],
             best.valid=ldply(.out, function(x) x$best.valid)[,2],
             best.pattern=ldply(.out, function(x) x$best.pattern)[,2])
}
