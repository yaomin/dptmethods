normalize.pmeans <-
function(count, countab, bg.correct=F, out.int=F) {
  ##browser()
  cnt.ex <- extend.tbs(countab)
  cnt.norm <- norm.tbs(cnt.ex, bg.correct=bg.correct)
  norm.sig <- as.data.frame(sapply(seq(ncol(count)), function(i) normalize.signal(count[,i],
                                                                                  cnt.ex[[i]],
                                                                                  cnt.norm[[i]],
                                                                                  out.int=out.int)))
  row.names(norm.sig) <- row.names(count)
  names(norm.sig) <- names(count)
  norm.sig
}
