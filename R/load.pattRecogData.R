load.pattRecogData <-
function(from, chr, pos=-1, pos.out=parent.frame()) {
  if(missing(chr)) chr <- get('chr', pos=pos)
  if(missing(from)) from <- get.ws.path('pattRecog', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','pattRecog')]],chr,"RData", sep=".")),
       envir=as.environment(pos.out))
}
