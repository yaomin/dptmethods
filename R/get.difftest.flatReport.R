#' @export
get.difftest.flatReport <-
function(path,
                                    chrs,
                                    winsize,
                                    pos=-1
                                    )
{
  ##browser()
  dpt.options <- get("dpt.options", pos=pos)
  ws <- 'diffTest'
  if(missing(path)) path <- get.ws.path(ws, pos=pos)
  if(missing(chrs)) chrs <- get.ws.chrs(ws, pos=pos)
  report.flat <- NULL
  for (chr in chrs) {
    cat("loading",dpt.options[[c("output","diffTest")]],"data:",chr)
    load.diffTestData(path, chr, pos=pos)
    cat("\t done\n")
    report.flat <- rbind(report.flat,
                         flattern.testReport(wins.test.report,
                                             chr,
                                             winsize))
  }
  report.flat
}
