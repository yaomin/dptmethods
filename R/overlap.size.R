overlap.size <-
function(R1,R2) {
  ## R1 and R2 are expected to be matrix with two columns and same number of rows.
  
  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  
  delta1 <- R1[,2]-R2[,1]
  delta2 <- R1[,1]-R2[,2]
  olp <- delta1*delta2<0
  deltap <- apply(cbind(delta1, delta2, R1[,2]-R1[,1], R2[,2]-R2[,1]), 1, 
                  function(x) min(abs(x)))
  deltap
}
