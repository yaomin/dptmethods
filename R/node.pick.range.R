node.pick.range <-
function(anode, id.start=1,tol=1e4, debug=F){
  rests <- NULL
  i <- id.start
  getmembers <- function(X) {
    if(!is.leaf(X)) {
      output <- F
      pr0 <- attr(X, 'Pr')
      pr1 <- attr(X[[1]], 'Pr')
      pr2 <- attr(X[[2]], 'Pr')
      lf1 <- is.leaf(X[[1]])
      lf2 <- is.leaf(X[[2]])
      if(debug) cat(pr0,pr1,pr2,lf1,lf2)
      if(attr(X, 'Sig')) {
        if(lf1 & !lf2 & pr2 > pr0) output <- T
        else if(!lf1 & lf2 & pr1 > pr0) output <- T
        else if(lf1 & lf2) output <- T
        else if(!lf1 & !lf2 & (pr0-min(pr1, pr2))<tol) output <- T
      }

      if(debug) cat('-->',output,'\n')

      if(output) {
        this.wins <- sort(labels(X))
        this.ID <- rep(i, length(this.wins))
        this.pr <- rep(attr(X, "Pr"), length(this.wins))
        this.res <- cbind(this.ID, this.wins, this.pr)
        rests <<- rbind(rests,this.res)
        i <<- i+1
        
      }
      else {
        for (jj in seq(length=length(X))) {
          X[[jj]] <- getmembers(X[[jj]]) 
        }
      }
    }
    X
  }
  getmembers(anode)
  if(!is.null(rests)) {
    rests <- as.data.frame(rests)
    names(rests) <- c('ID','Win','Pr')
  }
  rests
}
