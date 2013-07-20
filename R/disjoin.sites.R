#' @export
disjoin.sites <-
function(site.rl,
	 gaps.within.rl,
	 type=c("any", "start", "end", "two.ends","within","equal")) {
    #    browser()
    type <- match.arg(type)
    if(type =="two.ends") {
	gaps.mtched <- gaps.within.rl[!gaps.within.rl%within%site.rl]
  } else {
    gaps.mtched <- subsetByOverlaps(gaps.within.rl, site.rl, type=type)
  }
  ## site with gaps on either end
  site.rl.gaps <- site.rl[site.rl%over%gaps.mtched]
  ## site with nogaps on either end
  site.rl.nogaps <- site.rl[!site.rl%over%gaps.mtched] 
  ## disjoin gaps and sites with gaps
  .disjoin <- disjoin(combine.2rl(site.rl.gaps, gaps.mtched))
  ## remove the gaps in sites with end gaps 
  site.rl.gaps.cleaned <- .disjoin[!(.disjoin%over%gaps.mtched)]
  out <- combine.2rl(site.rl.nogaps, site.rl.gaps.cleaned)
  out
}
