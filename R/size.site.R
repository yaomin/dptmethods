size.site <-
function(sites) {
  ## browser()
  output <- NULL
  sites.unique.pattern <- unique(sites$pattern)
  for (i in seq(along=sites.unique.pattern)) {
    this.sites <- subset(sites, subset=pattern==sites.unique.pattern[i])
    sites.Win <- as.numeric(this.sites$Win)
    sites.ID <- as.character(this.sites$ID)
    ##sites.pattern <- sites$pattern
    this.output <- tapply(sites.Win,
                          sites.ID,
                          function(x) max(x)-min(x),
                          simplify=T)
    idx <- match(unique(sites.ID), names(this.output))
    output <- rbind(output, data.frame(ID=unique(sites.ID),
                                       size=this.output[idx],
                                       pattern=sites.unique.pattern[i],
                                       stringsAsFactors=F))
  }
  
  output
}
