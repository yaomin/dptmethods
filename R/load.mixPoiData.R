#' @export
load.mixPoiData <-
function(from, chr, pos=-1, pos.out=parent.frame()) {
  if(missing(chr)) chr <- get('chr', pos=pos)
  if(missing(from)) from <- get.ws.path('mixPoisson', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from,paste(dpt.options[[c('output','mixPoisson')]],chr,"RData", sep=".")),
       envir=as.environment(pos.out))
}
