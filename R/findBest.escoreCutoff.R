findBest.escoreCutoff <-
function(escores, cuts=seq(-2,2,0.01), cutoff=0.1) {
  n.tot <- nrow(escores)[1]
  miss.n <- sapply(cuts, function(k, dat) sum(rowSums(dat > k) > 1), dat=escores)
  miss.pct <- miss.n/n.tot
  cuts[which(miss.pct<=cutoff)[1]]
}
