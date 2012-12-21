est.beta.para <-
function(x) {
  ## Moment estimate of parameters
  u <- mean(x)
  v <- var(x)
  a <- (1-u)*u^2/v-u
  b <- (1-u)*(u*(1-u)/v-1)
  list(a=a, b=b)
}
