attach.pmeans.cutoff <-
function(sites.df,
                                 score.ecdf,
                                 by.chr=F,
                                 pattern.col="pattern",
                                 size.col="size",
                                 score.col="pmeans",
                                 chr.col="chr",
                                 winsize=50)
{
  ##browser()
  patts <- unique(sites.df[,pattern.col])
  chrs <- unique(sites.df[, chr.col])
  out.site <- NULL
  if(by.chr) {

  } else {
    for(ipatt in patts) {
      sites.df.patt <- subset(sites.df, subset=pattern==ipatt)
      bf.ecdf.patt <- score.ecdf[[1]]
      bfP <- apply(sites.df.patt, 1, function(x) {
        ##print(x);
        this.ecdf <- bf.ecdf.patt[[as.integer(as.numeric(x[size.col])/winsize)]];
        1-this.ecdf(as.numeric(x[score.col]))})
      iout <- cbind(sites.df.patt,bfP)
      names(iout) <- c(names(sites.df.patt), "score.P")
      out.site <- rbind(out.site, iout)
    }
  }
  out.site
}
