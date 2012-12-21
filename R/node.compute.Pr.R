node.compute.Pr <-
function(anode, n, p, k.cutoff=100,win.size=100) {
  ##browser()
  if(!is.leaf(anode)) {
    ewins <- as.numeric(labels(anode))/win.size
    k <- (max(ewins)-min(ewins))
    r <- attr(anode, 'members')
    ##cat(r,'\t',k,'\t',n,'\t',p,'\n')
    if(k>=r/p | k>k.cutoff) 1 else 1-Pr.S.up(r,k,n,p)}
  else NA
}
