load.init.ws <-
function(path, pos=-1) {
  if(missing(path)) {
    load(file.path(get.ws.path("root", pos=pos),".RData"))
  }
  else load(file.path(path, ".RData"))
}
