which.max.by <-
function(by, score) {
  
  .id <- seq(length(score))
  score.order <- order(score, na.last=T,
                       decreasing=T)
  .select <- vaggregate(.id[score.order],
                        factor(by[score.order]),
                        function(x) x[1])
  .select
}
