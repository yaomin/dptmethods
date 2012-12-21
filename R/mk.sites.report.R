mk.sites.report <-
function(sites,
                            index.col,
                            win.col="Win",
                            win.score="Win.score",
                            score.fun=c("max","mean"),
                            winsize=50,
                            ...)
{
  rangeFun <- function(x) {
    r.x <- range(as.integer(x[,win.col]))
    s.x <- sapply(score.fun, function(y) do.call(y, list(x[, win.score],...)))
    c(r.x[1], r.x[2]+winsize, diff(r.x)+winsize, s.x)
  }
  
  ret <- ddply(sites, index.col, rangeFun)
  s.id <- apply(ret[,index.col], 1, paste, collapse="-", sep="")
  ret <- cbind(ret, s.id)
  
  scorenames <- if(length(score.fun) >1) {
    c("score", paste("score", seq(length(score.fun)-1), sep="."))
    } else "score"
  names(ret) <- c(index.col, c("start","end","size", scorenames, "siteID"))
  ret
}
