get.sharedSource <-
function(url, pos=-1) {
  src <- getURL(url)
  ##ws.env <- if(pos==-1) globalenv() else as.environment(pos)
  eval(parse(text = src), as.environment(pos))
}
