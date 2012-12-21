merge.fimo.site <-
function(f, s) {
  
  s.1 <- subset(s, subset=hf.fimo)
  s.2 <- subset(s, subset=!hf.fimo)
  row.names(s.1) <- as.character(s.1$siteID)
  sf.1 <- s.1[as.character(f$V2),]
  full.1 <- cbind(sf.1,f)
  full.1.names <- names(full.1)

  sf.2 <- matrix(NA, nrow=nrow(s.2), ncol=ncol(f))
  full.2 <- cbind(s.2, sf.2)
  names(full.2) <- full.1.names

  full <- rbind(full.1, full.2)
  row.names(full)<- seq(nrow(full))
  full
}
