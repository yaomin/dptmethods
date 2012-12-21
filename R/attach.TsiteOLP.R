attach.TsiteOLP <-
function(sites.chr, Tsites) {
  
  stopifnot(length(unique(sites.chr$chr)) ==1)
  chrom <- as.character(unique(sites.chr$chr))
  cat(chrom, "\n")
  ##browser()
  overlap <- vector("logical", nrow(sites.chr))
  whichtsite <- overlapProp.tsite <- vector("numeric", nrow(sites.chr))
  overlapProp.site<- overlapSize <- vector("numeric", nrow(sites.chr))
  sites <- subset(sites.chr, select=c("start","end"))
  tsites <- subset(Tsites, subset=chr==chrom, select=c("start","end"))
  if(nrow(tsites) > 1) {
    olp <- any.overlap.sites(tsites, sites, e=0)
    ##olp$where <- as.integer(as.character(olp$where))
    olp.where <- as.numeric(unlist(strsplit(paste(as.character(olp[olp$overlap,"where"]),
                                                  collapse=",",
                                                  sep=""),
                                            split=","))
                            )

    olp.noNA <- subset(olp, subset=overlap)
    new.olp.noNA <- ldply(apply(cbind(tsite=row.names(olp.noNA),olp.noNA),
                                1,
                                function(x) data.frame(cbind(x[1],
                                                             as.numeric(unlist(strsplit(x[4], split=",")))
                                                             )
                                                       )
                                ),
                          function(y) y
                          )
    new.olp.noNA$X1 <- as.numeric(as.character(new.olp.noNA$X1))
    new.olp.noNA$X2 <- as.numeric(as.character(new.olp.noNA$X2))
    
    tsite.olp <- tsites[new.olp.noNA$X1,]
    ##sites.olp <- sites[olp[olp$overlap, "where"],]
    sites.olp <- sites[new.olp.noNA$X2,]
    ##browser()
    
    ##overlap[olp.noNA$where] <- T
    overlap[new.olp.noNA$X2] <- T
    ##overlap[olp.where] <- T
    
    ##whichtsite[olp.noNA$where] <- row.names(tsites)[olp$overlap]
    whichtsite[new.olp.noNA$X2] <- row.names(tsites)[new.olp.noNA$X1]
    olp.prop <- overlap.proportion(tsite.olp, sites.olp)
    ##overlapProp.tsite[olp.noNA$where] <- olp.prop$p1
    overlapProp.tsite[new.olp.noNA$X2] <- olp.prop$p1
    ##overlapProp.site[olp.noNA$where] <- olp.prop$p2
    overlapProp.site[new.olp.noNA$X2] <- olp.prop$p2
    overlapSize[new.olp.noNA$X2] <- overlap.size(tsite.olp, sites.olp)
  }
  out <- cbind(sites.chr,
               overlap=overlap,
               whichTsite=whichtsite,
               overlapProp.Tsite=overlapProp.tsite,
               overlapProp.site=overlapProp.site,
               overlapSize=overlapSize)
  out
}
