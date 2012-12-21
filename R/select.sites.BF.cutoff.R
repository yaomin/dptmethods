select.sites.BF.cutoff <-
function(sites.df,
                                   bf.cutoff,
                                   pattern.col="pattern",
                                   size.col="size",
                                   score.col="score",
                                   chr.col="chr",
                                   winsize=50
                                   )
{
  ## Assumption: sites.df
  ## > head(sites.df)
  ##     .id   chr pattern     ID    start      end size    score  score.1
  ## 1 chr10 chr10      P1  1.100 31360650 31361400  750 15.19609 12.31469
  ## 2 chr10 chr10      P1 1.1000 14603400 14603450   50 14.20799 14.20799
  ## 3 chr10 chr10      P1 1.1001 14625350 14625400   50 14.11777 14.11777
  ## 4 chr10 chr10      P1 1.1002 14646250 14646300   50 14.30350 14.30350
  ## 5 chr10 chr10      P1 1.1003 14681400 14681450   50 14.48476 14.48476
  ## 6 chr10 chr10      P1 1.1004 14730350 14730400   50 13.88597 13.88597
  ##            siteID
  ## 1  chr10-P1-1.100
  ## 2 chr10-P1-1.1000
  ## 3 chr10-P1-1.1001
  ## 4 chr10-P1-1.1002
  ## 5 chr10-P1-1.1003
  ## 6 chr10-P1-1.1004
  ##
  ##browser()
  patts <- unique(sites.df[,pattern.col])
  chrs <- unique(sites.df[, chr.col])
  out.sites <- NULL
  if("chr" %in% names(bf.cutoff)) {
    for (ichr in chrs) {
      for (ipatt in patts) {
        isites <- subset(sites.df, subset=sites.df[,pattern.col]==ipatt & sites.df[, chr.col]==ichr)
        ibfcutoff <- subset(bf.cutoff, subset=chr==ichr & pattern == ipatt)
        out.isites <- subset(isites, subset=isites[,score.col]>ibfcutoff[isites[,size.col]/50, 'cutoff'])
        out.sites <- rbind(out.sites, out.isites)
      }
    }
  } else {
    for (ipatt in patts) {
      isites <- subset(sites.df, subset=sites.df[,pattern.col]==ipatt)
      ibfcutoff <- subset(bf.cutoff, subset=pattern == ipatt)
      out.isites <- subset(isites, subset=isites[,score.col]>ibfcutoff[isites[,size.col]/50,"cutoff"])
      out.sites <- rbind(out.sites, out.isites)
    }
  }
  out.sites
}
