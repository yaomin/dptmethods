win.overlap.pair <-
function(wins.sel.pair) {
  wins.sel.1 <- wins.sel.pair[[1]]
  wins.sel.2 <- wins.sel.pair[[2]]
  win.int <- intersect(wins.sel.1$wins, wins.sel.2$wins)
  del.in.2 <- wins.sel.2[win.int,"probs"] >= wins.sel.1[win.int,"probs"]
}
