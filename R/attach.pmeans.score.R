attach.pmeans.score <-
function(sites.js.bfcut.1,
                                index.col,
                                sites.js.ex,
                                pmeans.norm,
                                sample.select) {
  res1 <- pmeans.norm[sites.js.ex$Win,sample.select, drop=F]
  res2 <- data.frame(sites.js.ex, res1)
  names(res2) <- c(names(sites.js.ex), "pmeans")
  res4 <- ddply(res2, index.col, function(x) mean(x$pmeans))$V1
  res5 <- cbind(sites.js.bfcut.1, res4)
  names(res5) <- c(names(sites.js.bfcut.1), "pmeans")
  res5
}
