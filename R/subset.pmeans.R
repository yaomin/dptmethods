subset.pmeans <-
function(pmeans, subset) {
  list(mean=pmeans$mean[,subset],
       var=pmeans$var[,subset])
}
