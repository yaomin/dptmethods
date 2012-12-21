extract.info.fromQuESTPeakOutput <-
function(peakNames) {
  process.name <- function(aname) {
    ##browser()
    nmfull <- unlist(strsplit(aname, split="~"))
    nmout.1 <- nmfull[c(1,2,3)]
    nmout.2 <- unlist(strsplit(nmfull[9], split="-"))
    nmout.3 <- nmfull[c(17,19)]
    nmout <- c(nmout.1, nmout.2, nmout.3)
  }
  fullInfo <- ldply(sfClusterApplyLB(peakNames, process.name), function(x) x)
  row.names(fullInfo) <- peakNames
  fullInfo <- cbind(fullInfo, peakNames)
  names(fullInfo) <- c("ID", "chr","center","start","end","qv","rank","siteID")
  class(fullInfo$center) <- class(fullInfo$start) <- class(fullInfo$end) <- "numeric"
  class(fullInfo$qv) <- class(fullInfo$rank) <- "numeric"
  fullInfo
}
