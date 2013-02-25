#' @export
rangedSites <-
function(sites, winsize) {
  
  sites$Win <- as.integer(sites$Win)  
  sites.win <- RangesList(lapply(split(sites, sites$pattern),
                                 function(x) IRanges(start=x$Win, width=winsize)))

  s.win <- sites$Win
  s.pattern <- sites$pattern
  s.id <- sites$ID
  s.min <- tapply(s.win, list(s.id,s.pattern), min)
  s.max <- tapply(s.win, list(s.id,s.pattern), max)
  s.ends <- sapply(seq(ncol(s.min)),
                   function(icol,x,y) {idx <- is.na(x[,icol]);
                                       cbind(x[!idx,icol], y[!idx,icol]) },
                   x=s.min,
                   y=s.max)

  sites.startend <- lapply(s.ends, function(x) cbind(x[,1],x[,2]+winsize-1))
  names(sites.startend) <- dimnames(s.min)[[2]]
  
  sites.site <- RangesList(lapply(sites.startend,
                                  function(x) IRanges(start=x[,1],end=x[,2])))
  list(sites.win=sites.win,
       sites.site=sites.site)
}
