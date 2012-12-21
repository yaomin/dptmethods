get.winsize.input <-
function(winsize) {
  wsize <- as.numeric(gsub("winsize","",
                           gsub('-','',grep('winsize[:digit:]*', commandArgs(), v=T))))
  if(length(wsize) == 0) wsize <- 100
  if(missing(winsize)) wsize else winsize
}
