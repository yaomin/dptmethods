compete.patt <-
function(sitelist) {
  ##browser()
  olp <- findOverlaps(sitelist[[1]], sitelist[[2]])
  if(length(olp) >1) stop("Requre RangesMatchingList of length 1 between two sites")
  mm <- as.matrix(olp)
  if(nrow(mm)>0) {
    ## valid: scores in [[1]] > [[2]]
    s1 <- score(sitelist[[1]])[mm[,1]]
    s2 <- score(sitelist[[2]])[mm[,2]]
    r1 <- s1>s2
    if(sum(r1)==length(r1)) {
      sitelist[[1]] <- sitelist[[1]]
      sitelist[[2]] <- sitelist[[2]][-(mm[,2]),]
    } else if(sum(!r1)==length(r1)) {
      sitelist[[1]] <- sitelist[[1]][-(mm[,1]),]
      sitelist[[2]] <- sitelist[[2]]
    } else {
      sitelist[[1]] <- sitelist[[1]][-(mm[,1][!r1]),]
      sitelist[[2]] <- sitelist[[2]][-(mm[,2][r1]),]
    }
  }
  sitelist
}
