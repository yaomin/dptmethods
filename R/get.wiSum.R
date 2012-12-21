get.wiSum <-
function(res) {
  ret <- lapply(res, function(x) if(!is.null(x)) colSums(x$w.win) else NULL)
  ret <-  t(as.data.frame(ret[!unlist(lapply(ret, is.null))]))
  row.names(ret) <- NULL
  ret
}
