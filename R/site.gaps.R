#' @export
site.gaps <- function(rs) {
  
  gaps.win <- with(rs,
                   gaps(sites.win,
                        start=min(start(sites.win)),
                        end=max(end(sites.win))))
  gaps.site <- with(rs, gaps(sites.site,
                             start=min(start(sites.site)),
                             end=max(end(sites.site))))
  
  gaps.within <- with(rs, intersect(sites.site, gaps.win))
  
  list(between=gaps.site,
       within=gaps.within)
}
