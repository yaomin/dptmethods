get.str.chr <-
function(chrStr) {
  chr <- gsub('-','',grep('chr', commandArgs(), v=T))
  if(missing(chrStr)) {
    if(length(chr) == 0) stop("chr string is missing.")
    else chr
  } else chrStr
}
