win.movingCount <-
function(ares, win.size=4) {
  pulse.n <- dim(ares)[1]
  win.n <- pulse.n - win.size + 1
  res <- sapply(seq(win.n),
                function(x, dat, win.size) sum(dat[x - 1 + seq(win.size),2]),
                dat = ares,
                win.size=win.size)
  cbind(V2=ares[seq(win.n),1], V4=res)
}
