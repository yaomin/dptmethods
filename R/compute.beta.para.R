compute.beta.para <-
function(res, d = c(1/3, 1/3, 1/3)) {
  wisum <- get.wiSum(res)
  d <- matrix(rep(d, each=nrow(wisum)), ncol= length(d))
  wd <- wisum+d
  cbind(wd[,3], wd[,1]+wd[,2])
}
