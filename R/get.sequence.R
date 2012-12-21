get.sequence <-
function(sites, ref, winsize,ext=500) {
  ### sites -  is the general sites format produced at different stages of
  ###  DPT analysis. It requires to have columns: chr, pattern, ID,
  ###  start, end.
  ### ref - reference genome.
  outseq.names <- apply(sites, 1, function(x) paste(x[c("chr", "pattern", "ID")],
                                               collapse="-", sep=""))
  outseq <- subseq(ref, start=sites[1,'start']-ext, end=sites[1,'end']+winsize-1+ext)
  for(i in seq(nrow(sites))[-1]) {
    outseq <- c(outseq, subseq(ref, start=sites[i,'start']-ext, end=sites[i,'end']+winsize-1+ext))
  }
  names(outseq) <- outseq.names
  outseq
}
