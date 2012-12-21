flattern.testReport <-
function(test.report,
                                chr,
                                winsize,
                                colnames)
{
  ## test.report can either be wins.test or wins.test.report
  ##browser()
  ### Assume the same number of columns in the contrasts of test.report
  if(missing(colnames)) {
    ##n.ctr <- ncol(attr(test.report[[1]][[1]], "contrast")) (Not as reliable)
    n.ctr <- nrow(attr(test.report[[1]][[1]], "contrast"))-1
    seq.ctr <- seq(0,n.ctr)
    colnames=c("start","end","ID",
      paste("value",seq.ctr, sep="."),
      paste("p.value",seq.ctr, sep="."))
  }
  
  output <- NULL
  for(i in seq(along=test.report)) {
    for(j in seq(along=test.report[[i]])) {
      output.ij <- cbind(chr,
                         pattern=code.pattern(attr(test.report[[i]][[j]],'pattern')),
                         contrast=code.contrast(attr(test.report[[i]][[j]],'contrast')),
                         test.report[[i]][[j]]
                         )
      names(output.ij) <- c("chr","pattern","contrast",colnames)
      output <- rbind(output, output.ij)
    }
  }
  ###: CREATED LAST WIN PROBLEM output[,'end'] <- output[,'end']+winsize
  ##test.report.names <- names(test.report[[i]][[j]]) ## purposely
  ##names(output) <- c('chr','pattern','contrast',test.report.names)
  output
}
