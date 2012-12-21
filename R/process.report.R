process.report <-
function(chr, winsize=100) {
  regionPatt.file <- paste(region.path,
                           paste(chr,"consolidated_wins.RData", sep="."),
                           sep="/")
  mixPoi.file <- paste(mixPoi.path,
                       paste("mixPoi",chr,"RData", sep="."),
                       sep="/")
  diff.file <- paste(diff.path,
                     paste("differential",chr,"RData", sep="."),
                     sep="/")
  load(regionPatt.file)
  load(mixPoi.file)
  load(diff.file)

  patts <- unique(sites.js$pattern)
  patt.n <- length(patts)

  cat("*** ",chr," ***","\n")
  sel.test.report <- NULL
  sel.read.report <- NULL
  for (i in seq(patt.n)[-7]) {
    for(j in seq(length(events.ctr[[i]]))) {
      sel.test.report.0 <- select.sites(all.test.report[[i]][[j]],1,1,0,0,win.size=winsize)
      if(nrow(sel.test.report.0) > 0) {
        sel.read.report.0 <- selReads.report(row.names(sel.test.report.0),
                                             subset(sites.js, subset=pattern==patts[i]),
                                             reads)
  
        sel.test.report.0 <- data.frame(chr=chr,
                                        pattern=as.character(patts)[i],
                                        contrast=paste("c",
                                          paste(events.ctr[[i]][[j]][,1],
                                                collapse=""),
                                          sep=""),
                                        ID=row.names(sel.test.report.0),
                                        sel.test.report.0)
        sel.read.report.0 <- data.frame(chr=chr,
                                        pattern=as.character(patts)[i],
                                        contrast=paste("c",
                                          paste(events.ctr[[i]][[j]][,1],
                                                collapse=""),
                                          sep=""),
                                        sel.read.report.0)
        
        sel.test.report <- rbind(sel.test.report, sel.test.report.0)
        
        sel.read.report <- rbind(sel.read.report, sel.read.report.0)
        cat(as.character(patts[i]),"\t",j,"\t", nrow(sel.test.report.0),"\n")
      }
    }
  }
  list(sites=sel.test.report, reads=sel.read.report)
}
