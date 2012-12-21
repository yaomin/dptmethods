get.sequence.chr <-
function(sites, ref, wisize,ext=500) {
  ### Assumption: sites and ref are refering to the same chromosome.
  n.site <- nrow(sites); stopifnot(n.site>0)
  outseq <- subseq(ref,
                   start=sites[1,'start']-ext,
                   end=sites[1,'end']+winsize-1+ext)
  for(i in seq(nrow(sites))[-1]) {
    outseq <- c(outseq, subseq(ref,
                               start=sites[i,'start']-ext,
                               end=sites[i,'end']+winsize-1+ext))
  }
  ## outseq <- apply(sites, 1, function(x) subseq(ref,
  ##                                              start=as.integer(x[4]),
  ##                                              end=as.integer(x[5])-1))
  names(outseq) <- sites$siteID
  outseq
}
