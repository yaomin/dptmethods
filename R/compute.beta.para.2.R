compute.beta.para.2 <-
function(res) {
  a <- get.w3Chr(res)
  cbind(a=a, b=1-a)
}
