match.pattern.win <-
function(wins.sel, sites) {
  ##browser()
  pattern.win <- rep("", nrow(sites))
  sites.Win <- sites$Win
  wins.sel.names <- names(wins.sel)
  for(i in seq(along=wins.sel)) {
    idx <- sites.Win %in% wins.sel[[i]]$wins
    ## idx.noNA <- na.omit(idx)
    ## if(length(idx.noNA) > 0) {
    ##   print(length(idx.noNA))
    pattern.win[idx] <- paste(pattern.win[idx], wins.sel.names[i], sep="")
    ## }
  }
  pattern.win
}
