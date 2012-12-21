get.siteSig <-
function(sigres, msites, which.pattern, which.site) {
  Patt.idx <- grep(which.pattern,msites$PatternID)
  Win.idx <- which(msites$SiteID==which.site)
  wins <- msites[intersect(Patt.idx, Win.idx),'Win']
  sigres[as.character(wins),]
}
