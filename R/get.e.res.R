get.e.res <-
function(start, end, chr) {
  get.ws("pattRecog", what=c("e.res"), chr=chr)
  poi <- as.integer(rownames(e.res))
  .out <- e.res[poi>=start&poi<=end,]
  .out[order(as.integer(rownames(.out))),]
}
