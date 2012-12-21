selReads.report <-
function(id, wins, reads) {
  if('poi' %in% names(reads)) {
    reads.poi <- reads$poi
    reads <- reads[,-1]
    row.names(reads) <- reads.poi
  }

  id <- as.character(id)
  idx <- unlist(sapply(id, function(x,y) which(x==y), y=wins$ID))
  win.sel <- wins[idx,c('ID', 'Win')]
  read.sel <- reads[as.character(win.sel$Win),]
  res <- cbind(win.sel, read.sel)
  res
}
