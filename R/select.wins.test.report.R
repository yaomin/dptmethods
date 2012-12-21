select.wins.test.report <-
function(xl, method="fdr",
                                        cutoff=0.05,
                                        pvalue.col.list,
                                        effect.col.list,
                                        sign.names=c("Hypo","Hyper"))
{
  ## pvalue.col can be a vector: different elements in xl might take different p value cols.
  ##browser()
  n.list <- length(xl)
  output <- vector("list", length=n.list)

  ## if(length(pvalue.col.list) ==1) pvalue.col.list <- rep(pvalue.col.list, n.list)
  ## else stopifnot(length(pvalue.col.list) == n.list)

  ## if(length(effect.col) ==1) effect.col.list <- rep(effect.col.list, n.list)
  ## else stopifnot(length(effect.col.list) == n.list)
    
  for (i in seq(along = xl)) {
    pvalue.cols.i <- pvalue.col.list[[i]]
    effect.cols.i <- effect.col.list[[i]]
    output[[i]] <- vector("list", length(pvalue.cols.i))
    for (j in seq(along = pvalue.cols.i)) {
      
      pvalue.col.ij <- pvalue.cols.i[[j]]
      effect.col.ij <- effect.cols.i[[j]]
      xl.names.ij <- names(xl[[c(i,j)]])
      worker.ij <- function(k) {
        cat(i,"\t",j, "\t", k,"\n")
        ##pvalue.colname.j <- xl.names.ij[pvalue.col.j[k]]
        ##effect.colname.j <- xl.names.ij[effect.col.j[k]]
        ret <- select.wins.test.element(xl[[c(i,j)]],
                                        method = method,
                                        cutoff = cutoff,
                                        pvalue.col = pvalue.col.ij[k],
                                        effect.col = effect.col.ij[k],
                                        sign.names = sign.names)
        ret
      }
      output[[c(i,j)]] <- sapply(seq(along=pvalue.col.ij), worker.ij, simplify=F)
    }
  }
  output
}
