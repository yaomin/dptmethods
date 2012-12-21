get.w3Chr <-
function(res) {
  ret <- lapply(res, function(x) if(!is.null(x)) x$paraCRM$p[3] else NULL)
  c(unlist(ret))
}
