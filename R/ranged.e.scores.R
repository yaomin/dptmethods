ranged.e.scores <-
function(e.res, winsize) {
  RangedData(ranges=IRanges(start=as.integer(rownames(e.res)), width=winsize),
             e.res, space=NULL)
}
