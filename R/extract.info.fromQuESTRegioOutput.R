extract.info.fromQuESTRegioOutput <-
function(peakNames) {
  ##browser()
  process.name <- function(aname) {
    ##browser()
    nmfull <- unlist(strsplit(aname, split="~"))
    nmout.1 <- nmfull[c(1,2)]
    nmout.2 <- unlist(strsplit(nmfull[3], split="-"))
    nmout.3 <- nmfull[c(9,23,25)]
    nmout <- c(nmout.1, nmout.2, nmout.3)
  }
  fullInfo <- ldply(sfClusterApplyLB(peakNames, process.name), function(x) x)
  row.names(fullInfo) <- fullInfo[,1]
  fullInfo <- cbind(fullInfo, fullInfo[,1])
  names(fullInfo) <- c("ID", "chr","start","end","center","qv","rank","siteID")
  class(fullInfo$center) <- class(fullInfo$start) <- class(fullInfo$end) <- class(fullInfo$rank)<- "integer"
  class(fullInfo$qv) <- "numeric"
  fullInfo
}
