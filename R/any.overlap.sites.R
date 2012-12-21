any.overlap.sites <-
function(R1, R2, e=500, self=F) {

  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  R1.ex <- R1
  R1.ex[,1] <- R1[,1]-e
  R1.ex[,2] <- R1[,2]+e
  overlap.count <- rep(0, nrow(R1))
  overlap <- vector("logical", nrow(R1))
  overlap.where <- rep(NULL, nrow(R1))

  for (i in seq(nrow(R1))) {
    r1 <- R1.ex[i,]
    if(self) {
      overlap.v <- ((R2[-i,1]-r1[2]) * (R2[-i,2]-r1[1])) <0
      overlap.count[i] <- sum(overlap.v)
      overlap[i] <- overlap.count[i]>0
      overlap.where[i] <- paste(which(overlap.v), collapse=",")
    }
    else {
      overlap.v <- ((R2[,1] - r1[2]) * (R2[,2] - r1[1])) < 0
      overlap.count[i] <- sum(overlap.v)
      overlap[i] <- overlap.count[i]>0
      overlap.where[i] <- paste(which(overlap.v), collapse=",")
    }
    
  }
  data.frame(overlap=overlap, count=overlap.count, where=overlap.where)
}
