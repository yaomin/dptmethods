rmultin <-
function(p){	
  u <- runif(1, 0, 1)
  
  if (u<=p[1]){
    z <- c(1, 0, 0)
  } else {
    if (u<=(p[1]+p[2])){
      z <- c(0,1,0)	
    } else {
      z <- c(0,0,1)
    }
  }
  z
}
