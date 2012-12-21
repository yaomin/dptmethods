compute.minB.para <-
function(res,
                              pattern=c(0,0,0,1,1,1,1,1,1),
                              v=3)
{
  if(v==1) ab.1 <- compute.beta.para(res)
  else if(v==2) ab.1 <- compute.beta.para.2(res)
  else  ab.1 <- compute.beta.para.3(res)
  e.1 <- get.w3Chr(res)
  if(length(pattern) != nrow(ab.1) | length(pattern) != length(e.1)) stop("dimension mismatch!")
  ab.0 <- ab.1[,rev(seq(ncol(ab.1)))]
  e.0 <- 1 - e.1
  
  ab <- ab.1
  e <- e.1

  ab[pattern == 1,] <- ab.1[pattern == 1,]
  ab[pattern == 0,] <- ab.0[pattern == 0,]

  e[pattern == 1] <- e.1[pattern == 1]
  e[pattern == 0] <- e.0[pattern == 0]

  list(ab=ab, e=e)
}
