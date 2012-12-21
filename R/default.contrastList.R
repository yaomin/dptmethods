default.contrastList <-
function(patts) {
  ##browser()
  inpatt <- as.matrix(patts)
  out <- vector("list",ncol(inpatt))
  for(i in seq(dim(inpatt)[2])) {
    out[[i]] <- default.pattContrast(inpatt[,i])
  }
  unique(out)
}
