gen.poi <-
function(poi.data, sizes, win.size=100)
{ ## generate the positions
  ##browser()
  pois.start <- sample(poi.data,length(sizes))
  pois.end <- pois.start+sizes
  R1 <- cbind(pois.start, pois.end)
  Rex <- R1
  Rex[,1] <- R1[,1]-2000
  Rex[,2] <- R1[,2]+2000
  ##overlaps <- any.overlaps(R1, R1, e=2000)
  overlaps <- vector("logical", nrow(R1))
  for(i in seq(nrow(Rex))) {
    r1 <- Rex[i,]
    overlaps[i] <- any((Rex[-i,1]-r1[2]) * (Rex[-i,2]-r1[1]) <0)
  }
  R1 <- R1[!overlaps,]
  apply(R1, 1, function(x, by) seq.int(from=x[1], to=x[2], by=by), by=win.size)
}
