init.ws.dirs <-
function(what=c("log","src","tmp"), pos=-1) {
  create.ws.dir("root")
  ws.root.path <- get.ws.path('root', pos=pos)
  current.dir <- getwd()
  setwd(ws.root.path)
  for (ipath in what) {
    create.ws.dir(ipath, pos=pos)
  }
  setwd(current.dir)
}
