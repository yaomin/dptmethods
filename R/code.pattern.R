code.pattern <-
function(patt.mx) {
  if(is.null(dim(patt.mx))) patt.mx <- as.matrix(patt.mx)
  cpatt <- apply(patt.mx, 2, paste, collapse="")
  gsub("NA","*", cpatt)
}
