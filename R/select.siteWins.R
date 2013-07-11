select.siteWins <-
function(patt, wins, sites, chr) {
  .wins <- wins[wins%over%sites[[patt]]]
  .fac <- match(.wins, sites[[patt]])
  .dt <- data.frame(chr, patt, cbind(.fac, start(.wins)))
  names(.dt) <- c("chr", "pattern", "ID", "Win")
  .dt
}
