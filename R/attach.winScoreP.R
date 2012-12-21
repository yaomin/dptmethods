attach.winScoreP <-
function(sites, wins.sel) {
  out <- NULL
  for ( ipatt in names(wins.sel)) {
    patt.sub <- subset(sites, subset=pattern==ipatt)
    patt.score.p <- wins.sel[[ipatt]][patt.sub$Win,]
    out <- rbind(out, patt.score.p)
  }
  outall <- cbind(sites, out[,c('scores','probs')])
  names(outall) <- c(names(sites), c('Win.score', 'Win.P'))
  outall
}
