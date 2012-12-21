any.overlap.sites.allchr <-
function(R1,R2,
                                     e=500, self=F,
                                     chr.col="chr",
                                     start.col="start",
                                     end.col="end") {
  chrs <- unique(R1[,chr.col])
  overlapped.R1.sites <- NULL
  for (chr.i in chrs) {
    R1.sub <- subset(R1,
                     subset=R1[,chr.col]==chr.i,
                     select=c(start.col, end.col))
    R2.sub <- subset(R2, subset=R2[,chr.col]==chr.i,
                     select=c(start.col, end.col))
    olp <- any.overlap.sites(R1.sub,
                             R2.sub,
                             e=e,
                             self=self)
    row.names(olp) <- row.names(R1.sub)
    overlapped.R1.sites <- rbind(overlapped.R1.sites, olp)
    
  }
  overlapped.R1.sites
}
