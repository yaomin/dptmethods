clean.dupcontr.sites <-
function(sites,
                                 p.value.col,
                                 n.para=10)
{
  ## this seems to be the best I can do in R
  ##browser()
  siteIDs <- siteUID(sites)
  cat("# of site IDs:", length(siteIDs),"\n")
  cat("# of unique site IDs:", length(unique(siteIDs)),"\n")
  dupIdx <- duplicated(siteIDs)
  #dupIDs <- unique(siteIDs[table(siteIDs)>1])
  dupIDs <- unique(siteIDs[dupIdx])
  cat("# of unique site IDs with duplicates:", length(dupIDs), "\n")
  idx.dup <- siteIDs%in%dupIDs
  cat("# of sites with duplicates", sum(idx.dup), "\n")
  out.0 <- data.frame(.id=siteIDs[!idx.dup],sites[!idx.dup,])
  cat("# of sites with unique IDs", nrow(out.0), "\n")

  sites.dup <- sites[idx.dup,]
  sites.dup.siteIDs <- siteIDs[idx.dup]

  sf.fun.1 <- function(idx.sel) {
    sites.sel <- sites.dup[idx.sel,]
    out.sub <- by(sites.sel,
                  list(chr.ID=paste(sites.sel$chr,
                         sites.sel$pattern,
                         sites.sel$ID,                     
                         sep="-")),
                  function(x) {
##                    nr <- dim(x)[1]
##                    if(nr<2) x
##                    else x[which.min(x[,p.value.col]),]
                    x[which.min(x[,p.value.col]),]
                  })
    ldply(out.sub, function(x) x)  
  }
 
  sites.dup.cut <- cut(seq(length(dupIDs)), n.para)
  sites.dup.cut.levels <- levels(sites.dup.cut)

  sf.fun.2 <- function(cutlevel) {
    this.idx <- which(sites.dup.cut == cutlevel)
    next.idx <- which(sites.dup.siteIDs%in%dupIDs[this.idx])
    sf.fun.1(next.idx)
  }

  out.2.list <- sfLapply(sites.dup.cut.levels, sf.fun.2)
  out.2 <- ldply(out.2.list, function(x) x)
  out <- rbind(out.0,out.2)
  cat("# of sites returned:",nrow(out),"\n")
  out
}
