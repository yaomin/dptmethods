match.patt.pattmx <-
function(patt, pattmx) {
  match(patt, paste("P", apply(pattmx,2,paste, collapse=""), sep=""))
}
