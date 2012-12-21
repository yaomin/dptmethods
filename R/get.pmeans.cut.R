get.pmeans.cut <-
function(pmeans.df,lambda.vec, cut.p=0.05, method=2) {
  stopifnot(nrow(pmeans.df)!=length(lambda.vec))
  res <- NULL
  for ( i in seq(lambda.vec)) {
    if(method!=2) {
      res <- rbind(res, get.postMean.cutoffs.aux(pmeans.df[,i], lambda.vec[i], cut.p=cut.p))
    }
    else {
      this.res <- get.postMean.cutoffs.aux2(pmeans.df[,i], lambda.vec[i], cut.p=cut.p)
      res <- rbind(res, c(this.res, this.res))
    }
  }
  res <- as.data.frame(res)
  names(res) <- c('lo', 'hi')
  res
}
