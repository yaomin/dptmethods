tabsps <-
function(str) {
  paste(str, paste(rep("\t", (nchar(str)%/%6<1)+1), collapse="", sep=""), sep="")
}
