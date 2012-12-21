write.sites.seqs <-
function(siteseqs, filepath) {
  stopifnot(is(siteseqs, 'list'))
  chrs <- names(siteseqs)
  if(length(siteseqs) >0) {
    cat("\n", chrs[1])
    write.XStringSet(siteseqs[[1]], filepath)
    cat("\tDONE\n")
  }
  if(length(siteseqs)>1) {
    for (i in seq(along = siteseqs)[-1]) {
      cat(chrs[i])
      write.XStringSet(siteseqs[[i]], filepath, append=T)
      cat("\tDONE\n")
    }
  }
}
