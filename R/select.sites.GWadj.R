select.sites.GWadj <-
function(sites.data,
                               reads.data,
                               method=c("qvalue","BY","bonferroni"),
                               cutoff=0.05,
                               pvalue.colname="pvalue",
                               sites.out.file="allsites",
                               reads.out.file="allreads",
                               out.path="Reports")
{
  method <- match.arg(method)
  out.sites.sel <- NULL
  out.reads.sel <- NULL
  if(method == "qvalue") {
    q.sites <- compute.qvalues(sites.data[,pvalue.colname],lambda=0)
    out.sites.adj <- cbind(sites.data, qvalue=q.sites)
    out.sites.adj.sel <- subset(out.sites.adj, subset=qvalue<cutoff)
  }
  else if(method =="BY" | method =="bonferroni") {
    adj.sites <- p.adjust(sites.data[,pvalue.colname], method=method)
    out.sites.adj <- cbind(out.sites, p.adj=adj.sites)
    out.sites.adj.sel <- subset(out.sites.adj, subset=p.adj<cutoff)
  }
  chrs <- unique(out.sites.adj.sel$chr)
  for (achr in chrs) {
    cat("*** ", as.character(achr), " ***", "\n")
    sites.sub <- subset(out.sites.adj.sel, subset=chr==achr)
    patts.sub <- unique(sites.sub$pattern)
    for (apatt in patts.sub) {
      cat(as.character(apatt),"\t")
      sites.sub.patt <- subset(sites.sub, subset=pattern==apatt)
      ctrs.sub <- unique(sites.sub.patt$contrast)
      for (actr in ctrs.sub) {
        cat(as.character(actr),"\t")
        ids <- subset(sites.sub.patt, subset=contrast==actr, select="ID", drop=T)
        cat(length(ids), "\t")
        out.reads.i <- subset(reads.data, subset=(chr==achr &
                                                  pattern == apatt &
                                                  contrast == actr &
                                                  ID %in% ids))
        out.reads.sel <- rbind(out.reads.sel, out.reads.i)
        cat(nrow(out.reads.i), "\n")
      }
    }
    
  }
  if(!is.null(sites.out.file)) {
    write.csv(out.sites.adj.sel, file=paste(out.path,
                                   paste(sites.out.file,"adj",method,cutoff,"csv",sep="."),
                                   sep="/"))
  }
  if(!is.null(reads.out.file)) {
    write.csv(out.reads.sel, file=paste(out.path,
                                 paste(reads.out.file,"adj",method,cutoff,"csv",sep="."),
                                 sep="/"))  
  }
  list(sites=out.sites.adj.sel,
       reads=out.reads.sel)
}
