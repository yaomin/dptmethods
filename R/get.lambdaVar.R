get.lambdaVar <-
function(x) {
  if(is.null(x)) res <- NULL
  else res <- var(x$chains$lambda.mcmc)
  res
}
