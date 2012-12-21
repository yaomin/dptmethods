compute.cutoff <-
function(r,n,p, cutoff=0.05) {
  k.seq <- seq(r,floor(r/p))
  k.idx <- which((1-Pr.S.up(r, k.seq, n, p) < cutoff))
  res <- k.seq[rev(k.idx)[1]]
  if(length(res)>0) res else NA
}
