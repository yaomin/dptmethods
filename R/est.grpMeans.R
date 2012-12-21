est.grpMeans <-
function(all.test.report) {
  out.list <- all.test.report
  for (i in seq(all.test.report)) {
    for (j in seq(all.test.report[[i]])) {
      this.df <- all.test.report[[i]][[j]]
      a <- this.df[,'beffect']
      b <- this.df[,'effect']
      c <- this.df[,'eqffect']
      this.ex <- data.frame(baseline=a - b - c,
                            grp1=a + 2 * c,
                            grp2=a + b - c)
      out.list[[i]][[j]] <- cbind(this.df, this.ex)
     }
  }
  out.list
}
