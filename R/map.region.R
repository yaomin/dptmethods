#' @export
map.region <-
function(test.report, gap) {
  ##browser()
### test.report has chr, start, and end columns
### gap is the max gap that will get joined, so the minimum gap after joining will be gap+1
  names(test.report) <- replace.names(names(test.report), "chr", "space")
  .rd <- as(test.report, "RangedData")
  .rl <- ranges(.rd)
  region.rl <- reduce(.rl, drop.empty.ranges=T, min.gapwidth=gap+1)
  ovlp <- as.matrix(findOverlaps(.rl, region.rl))
  r.id <- paste(.rd$space, "r",ovlp[,2], sep="_")
  .rd$r.id <- r.id
  out <- as.data.frame(.rd)
  names(out) <- replace.names(names(out), "space", "chr")
  out

}
