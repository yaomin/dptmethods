#' @export
format.sites4diffTest <-
function(wins, sites, chr) {
  .list <- lapply(names(sites[sapply(sites, length)>0]),
                  select.siteWins,
                  wins=wins,
                  sites=sites[sapply(sites, length)>0],
                  chr=chr)
  names(.list) <- names(sites[sapply(sites, length)>0])
  .out <- ldply(.list, function(x) x)[,-1]
  .out
}
