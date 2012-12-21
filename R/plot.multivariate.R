plot.multivariate <-
function(dist.cor, which=c(1,2), clust.method="complete",...)
{ 
  ## 1 for cluster and 2 for mds
  if(1 %in% which) {
    clust <- hclust(dist.cor, method=clust.method)
    plot(as.dendrogram(clust),...)
    
  }
  if(2 %in% which) {
    mds <- cmdscale(dist.cor,k=4)
    plot(mds[,1], mds[,2], bty="n", xlab="Dimension 1", ylab="Dimension 2",...)
    text(mds[,1], mds[,2], labels=row.names(mds), pos=1, cex=0.7)
    plot(mds[,1], mds[,3], bty="n", xlab="Dimension 1", ylab="Dimension 3",...)
    text(mds[,1], mds[,3], labels=row.names(mds), pos=1, cex=0.7)
    plot(mds[,2], mds[,3], bty="n", xlab="Dimension 2", ylab="Dimension 3",...)
    text(mds[,2], mds[,3], labels=row.names(mds), pos=1, cex=0.7)
  }
}
