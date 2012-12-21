map.regs.chr <-
function(trf.sig,
                         trf,
                         chrom,
                         methylType,
                         patts=c("**1","*1*","1**"),
                         cutoff=5,
                         winsize=100)
{
  
  trf.sig.chr <- subset(trf.sig, subset=chr==chrom & MethylType==methylType)
  trf.sig.ord <- trf.sig.chr[order(trf.sig.chr$start),]

  out <- trf.sig.ord[,c('ID','pattern')]
  for(i in seq(patts)) {
    patt <- patts[i]
    trf.s.part <- subset(trf,
                         subset= (pattern == patt &
                                  chr==chrom)
                         )
    trf.s <- cbind(trf.s.part, size=with(trf.s.part, end - start))
    trf.s.c <- subset(trf.s, subset = size > cutoff*winsize)
    trf.s.c.ord <- trf.s.c[order(trf.s.c$start),]

    idx <- any.overlap.sites(trf.sig.ord[,c('start','end')],
                             trf.s.c.ord[,c('start','end')],
                             e = 0,
                             self=F)
    unique.overlap <- idx$count==1
    trf.sig.idx <- seq(nrow(idx))[unique.overlap]
    ## trf.idx <- as.numeric(unlist(strsplit(paste(idx[idx$overlap,]$where,
    ##                                             collapse=","),
    ##                                       split=",")))
    trf.idx <- as.numeric(as.character(idx[unique.overlap,'where']))
    
    out.here <- matrix(rep(NA, 2*nrow(trf.sig.ord)), ncol=2)
    out.here[trf.sig.idx,1] <- trf.s.c.ord[trf.idx, 'ID']
    out.here[trf.sig.idx,2] <- as.character(trf.s.c.ord[trf.idx, 'pattern'])
    
    out <- cbind(out, out.here)
  }
  names(out) <- c('s.ID', "s.Pattern",
                  sapply(seq(patts), function(x) paste(c('r.ID','r.Pattern'),x, sep=".")))
  out
}
