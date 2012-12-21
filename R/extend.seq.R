extend.seq <-
function(seqdat, start=NULL, end=NULL, by=100) {
  
  if(is.null(start)) if(nrow(seqdat) > 0) NA else start=min(seqdat[,1])
  if(is.null(end)) if(nrow(seqdat) > 0) NA else end = max(seqdat[,1])
  V2 = seq(from=start, to=end, by = by)
  V4 = rep(0, length(V2))
  exted <- as.data.frame(cbind(V2,
                               V4))

  if(nrow(seqdat) > 0) exted[match(seqdat[,1],exted$V2),2] <- seqdat[,2]
  exted
}
