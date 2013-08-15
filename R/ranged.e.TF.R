#' @export
ranged.e.TF <-
function(e.res, winsize, w.sel) {
    e.res[] <- F
    e.wins <- row.names(e.res)
    for (patti in names(e.res)){
	sel.i <- w.sel[[patti]]$wins
	e.res[[patti]][e.wins%in%sel.i] <- T
    }
    RangedData(ranges=IRanges(start=as.integer(rownames(e.res)), width=winsize),
	       e.res)
}
