win.den <-
function(wins.sel.i, win.size=100) {

  #sel.wins <- as.numeric(wins[e>e.cut])
  #sel.wins.df <- data.frame(sel.wins, row.names=sel.wins)
  sel.wins.df <- wins.sel.i[,2,drop=F]
  den.res <- as.dendrogram(hclust(dist(sel.wins.df/win.size), method='complete'))
  den.res.reorder <- reorder(den.res, sel.wins.df$wins, mean)
  den.res.reorder
}
