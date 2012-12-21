get.sharedSrc <-
function(str, ver="current", pos=-1) {
  pathstr <- paste("http://dl.dropbox.com/u/10421230/dpt", ver, sep="/")
  urlStr <- paste(pathstr, str, sep="")
  get.sharedSource(urlStr, pos=pos)
}
