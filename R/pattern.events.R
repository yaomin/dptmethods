pattern.events <-
function(wins.sel, events.comp) {
  ## events.comp
  # NA NA
  # NA 1
  # 1  1
  ## use NA for elements that are free of change
  ##browser()
  events.comp <- as.matrix(events.comp)
  events <- attr(wins.sel, "event")
  comp.list <- NULL
  out.list <- NULL
  for (i in seq(ncol(events.comp))) {
    rows2sel <- !is.na(events.comp[,i])
    events.i <- as.matrix(events[rows2sel,,drop=F])
    events.sel <- which(!is.na(match.mx(t(events.i),
                                        t(events.comp[rows2sel,i,drop=F])
                                        )))
    comp.list <- c(comp.list, list(events.sel))
  }
  for (j in seq(comp.list)) {
    this.comp <- NULL
    for(l in comp.list[[j]]) {
      this.comp <- rbind(this.comp, wins.sel[[l]])
    }
    out.list <- c(out.list, list(this.comp))
  }
  out.list.names <- paste("P",apply(events.comp, 2, paste, collapse=""), sep="")
  out.list.names <- gsub("NA", "-", out.list.names)
  names(out.list) <- out.list.names
  ### take care the null list
  nonullidx <- !sapply(out.list, is.null)
  out.list <- out.list[nonullidx]
  ###
  attr(out.list, 'event') <- events.comp[,nonullidx,drop=F]
  attr(out.list, 'nonull') <- nonullidx
  out.list
}
