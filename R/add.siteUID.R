add.siteUID <-
function(tst) {
  sid <- paste(tst$chr, tst$pattern, tst$ID, sep="-")
  out <- transform(s.id=sid, tst)
  row.names(out) <- sid
  out
}
