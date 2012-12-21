get.winsize <-
function(path, pos=-1) {
  if(missing(path)) winsize.path <- get.ws.path("joinSample", pos=pos)
  load(paste(winsize.path, "winsize.RData",sep="/"))
  winsize
}
