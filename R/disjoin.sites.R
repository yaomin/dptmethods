#' @export
disjoin.sites <-
function(site.rl,
                          gaps.within.rl,
                          type=c("any", "start", "end", "two.ends","within","equal")) {
  type <- match.arg(type)
  if(type =="two.ends") {
    mtch1 <- subsetByOverlaps(gaps.within.rl, site.rl, type="start")
    mtch2 <- subsetByOverlaps(gaps.within.rl, site.rl, type="end")
    gaps.mtched <- combine.2rl(mtch1, mtch2)
  } else {
    gaps.mtched <- subsetByOverlaps(gaps.within.rl, site.rl, type=type)
  }
  .disjoin <- disjoin(combine.2rl(site.rl, gaps.mtched))
  out <- setdiff(.disjoin, gaps.mtched)
  out
}
