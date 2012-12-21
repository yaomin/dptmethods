escore.comp <-
function(x, group, which.p) {
  min(tapply(x,
             group,
             function(x, which.p) quantile(x, p=which.p),
             which.p=which.p))
}
