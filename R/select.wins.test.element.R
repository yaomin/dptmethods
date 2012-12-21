select.wins.test.element <-
function(x, method="fdr",
                                        cutoff=0.05,
                                        pvalue.col,
                                        effect.col,
                                        sign.names=c("Hypo","Hyper"))
{
  ##browser()
  pvalue.colname <- names(x)[pvalue.col]
  effect.colname <- names(x)[effect.col]
  
  sel <- subset(x, subset=p.adjust(x[,pvalue.colname], 'fdr') < 0.05)
  sign.idx <- (sel[,effect.colname] > 0) + 1
  output <- data.frame(sel, sign=sign.names[sign.idx])[,
                              c('start','end','ID',
                                "X.Intercept..Value",
                                effect.colname,
                                pvalue.colname,
                                "sign")]
  names(output) <- c('start','end','ID',"baseline","effect", "p.value", "sign")
  output
}
