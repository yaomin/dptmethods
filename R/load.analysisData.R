load.analysisData <-
function(from, pos=-1,pos.out=parent.frame()) {
  if(missing(from)) from <- get.ws.path('report', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','analysis')]],"RData", sep=".")),
       envir=as.environment(pos.out))
}
