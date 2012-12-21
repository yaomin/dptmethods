remove.dupcontr <-
function(flat.test.report.afterqvalue,
                            pvalue.var="p.value.1") {
  ##browser()
### flat.test.report.afterqvalue is a RangedData class or data.frame
  .df <- as.data.frame(flat.test.report.afterqvalue[!is.na(flat.test.report.afterqvalue[[pvalue.var]]),])
  .df$seq <- seq(nrow(.df))
  .df$id <- with(.df, paste(chr, pattern, ID, sep="_"))
  seq.select <-  which.min.by(.df$id, .df[[pvalue.var]])
  .df[seq.select,]
  
}
