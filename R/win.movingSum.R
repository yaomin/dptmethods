win.movingSum <-
function(ares, win.size=50, by=50) {
  pulse.n <- length(ares)
  win.n <- pulse.n - win.size + 1
  win.seq <- seq(1,win.n,by=by)
  res <- sapply(win.seq,
                function(x, dat, win.size) sum(dat[x - 1 + seq(win.size)]),
                dat = ares,
                win.size=win.size)
  res.sum <- cumsum(res)
  data.frame(winseq=win.seq+win.size-1, n=res.sum, pct=res.sum/(win.seq+win.size-1))
}
