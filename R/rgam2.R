rgam2 <-
function(shape, rate){
  
  n <- length(shape)
  if (min(shape)<1) stop("Invalid shape in rgam().")
  
  d <- shape-1/3  
  c <- 1/(3*sqrt(d))
  
  x <- rnorm(n, 0, 1)
  u <- log(runif(n, 0, 1))
  v <- (1+c*x)*(1+c*x)*(1+c*x)
  
  ###YXU: accept <- ifelse(v>0, 1, 0)
  accept <- (v>0)+0
  tmp <- v>0
  cc <- 0.5*x[tmp]*x[tmp]+d[tmp]-d[tmp]*v[tmp] + d[tmp]*log(v[tmp])
  ###YXU: accept[tmp] <- ifelse(u[tmp]<=cc, 1, 0)
  accept[tmp] <- (u[tmp]<=cc)+0
  
  yc <- d*v
  y <- yc
  y[accept!=1] <- NA
  ##browser()
  ##y <- ifelse(accept==1, yc, NA)
  y <- y/rate
  ##y <- ifelse(y>lam, y, NA)
  y
}
