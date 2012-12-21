find.cutoff.objective <-
function(x, paras, cutoff) {
  (compute.minB.p(x=x, paras=paras)-cutoff)^2
}
