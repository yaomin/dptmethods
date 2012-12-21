#'  @export
construct.contrast.objs <-
function(events.wins, events.ctr) {
  unique.ctrs <- unique(unlist(lapply(events.ctr,
                                      function(x) lapply(x,
                                                         function(x) paste(apply(as.matrix(x),
                                                                                 2,
                                                                                 paste,
                                                                                 collapse=""),
                                                                           collapse="_")))))
  patts.str <- apply(events.wins, 2, paste, collapse="")
  sapply(unique.ctrs, function(x) list(contrast=x,
                                       pattern=patts.str,
                                       pvalue.var="p.value.1",
                                       effect.var="value.1"),
         simplify=F)
}
