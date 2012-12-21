#' @export
ranged.e.TF <-
function(e.res, winsize, cutoff=-1) {
 RangedData(ranges=IRanges(start=as.integer(rownames(e.res)), width=winsize),
            e.res > cutoff)
}
