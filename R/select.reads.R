#' @export
select.reads <-
function(reads.df,
                         sample.select) {  
  poi <- reads.df$poi
  s.id <- reads.df$s.id
  cols.rmv <- match(c("poi", "s.id"), names(reads.df))
  if(missing(sample.select)) {
    reads.sub <- reads.df[,-cols.rmv]
  } else {
    reads.sub <- reads.df[,-cols.rmv][,rep(sample.select,2)]
  }
  cbind(poi=poi, reads.sub, s.id=s.id)
}
