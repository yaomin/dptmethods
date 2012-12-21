code.contrast <-
function(contr.mx) {
  if(is.null(dim(contr.mx))) contr.mx <- as.matrix(contr.mx)
  colcontr <- apply(contr.mx, 2, paste, collapse="")
  paste(colcontr, collapse="_")
}
