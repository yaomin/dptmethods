win.search <-
function(atree, cut=4e5, win.size=100) {
  atree.l <- cut(atree,cut)$lower
  res.wins <- NULL
  id.count <- 1
  allWins <- labels(atree)
  for (i in seq(atree.l)) {
    np <- node.get.p.n(atree.l[[i]],win.size)
    ##if (i >1 & nrow(res.wins) >0) id.count <- max(as.numeric(res.wins$ID))+1
    acut <- node.pick.range(node.set.significance(node.set.Pr(atree.l[[i]],
                                                              np$n, np$p,
                                                              win.size=win.size)),
                            id.start=(if(i==1) id.count else max(as.numeric(res.wins$ID))+1))
                            ##id.start = id.count)
    res.wins <- rbind(res.wins, acut)
  }
  res.wins$Win <- as.character(res.wins$Win)
  res.wins$ID <- as.integer(res.wins$ID)
  res.wins$Pr <- as.numeric(as.character(res.wins$Pr))
  singleWins <- setdiff(allWins, res.wins$Win)
  singleWins.res <- data.frame(ID=seq(singleWins)+length(unique(res.wins$ID)),
                               Win=as.character(singleWins), Pr=rep(0, length(singleWins)),
                               stringsAsFactors=F)
  rbind(res.wins, singleWins.res)
}
