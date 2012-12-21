make.reportContrast <-
function(report.flat,
                                pvalue.var="p.value.1",
                                effect.var="value.1",
                                type=c(1,2))
{
  if(type ==1) {
    contrDF <- expand.grid(dimnames(with(report.flat, table(contrast, pattern))))
    contrList <- vector("list", nrow(contrDF))
    for(i in seq(length(contrList))) {
      contrList[[i]]$contrast <- as.character(contrDF$contrast[i])
      contrList[[i]]$pattern <- as.character(contrDF$pattern[i])
      contrList[[i]]$pvalue.var <- pvalue.var
      contrList[[i]]$effect.var <- effect.var
    }
  }
  else if(type==2) {
    acontr <- dimnames(with(report.flat, table(contrast, pattern)))
    acontr$pvalue.var <- pvalue.var
    acontr$effect.var <- effect.var
    contrList <- list(acontr)
  }
  contrList
}
