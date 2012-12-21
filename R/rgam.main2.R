rgam.main2 <-
function(shape, rate){
  y <- rgam2(shape, rate)  
  ## Fill in NAs
  id <- is.na(y)
  this.shape <- shape[id]
  this.rate <- rate[id]
  y[id] <- rgam2(this.shape, this.rate)		
  
  ## Fill in NAs
  id <- is.na(y)
  this.shape <- matrix(shape[id], ncol=1)
  this.rate <- rate[id]

  y[id] <- rgamma(length(this.shape), shape=this.shape, rate=this.rate)
  ## for (i in 1:length(id)){
  ##   this.i <- id[i]
  ##   tmp <- rgamma(1, shape=this.shape[i], rate=this.rate[i])
  ##   y[this.i] <- tmp
  ## }	
  y
}
