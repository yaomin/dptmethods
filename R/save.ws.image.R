save.ws.image <-
function(ws.root.path, pos=-1) {
  if(missing(ws.root.path)) ws.root.path <- get.ws.path('root',pos=pos)
  current.dir <- getwd()
  setwd(ws.root.path)
  save.image()
  setwd(current.dir)
}
