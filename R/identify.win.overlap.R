identify.win.overlap <-
function(wins.sel) {
  
  ow <- NULL
  seq.n <- seq(length(wins.sel))
  for (i in seq.n ) {
    ow <- c(ow, wins.sel[[i]]$wins)
  }

  o.wins <-  names(which(table(ow)>1))

  ow.p <- NULL
  for (i in seq.n) {
    ow.p <- cbind(ow.p, wins.sel[[i]][o.wins,"probs"])
  }

  zero.col <- which(colSums(attr(wins.sel, 'event')) == 0)
  if(length(zero.col) ==1) {
    ow.p.this <- ow.p[,-zero.col,drop=F]
    if(ncol(ow.p.this)>1) {
      ids <- apply(ow.p.this,1, which.min)
    } else {
      ids <- NULL
    }
  }
  else if (length(zero.col) ==0) {
    if(ncol(ow.p)>1) {
      ids <- apply(ow.p,1, which.min)
    } else {
      ids <- NULL
    }
  }
  else stop("Wrong number of columns for zero patterns!")
  ##ids.full <- apply(ow.p,1, which.min)
  ret <- vector("list", length(wins.sel))
  for (i in seq.n) {
    o.wins.i <- o.wins[ids != i]
    wins.sel.i <- wins.sel[[i]]
    ret[[i]] <- wins.sel.i[as.character(setdiff(wins.sel.i$wins, o.wins.i)),]
  }
  attributes(ret) <- attributes(wins.sel)
  ret
}
