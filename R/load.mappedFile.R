load.mappedFile <-
function(filepath, cat.pre="", cat.after="\t") {
  cat(cat.pre, basename(filepath),cat.after)
  x <- scan(filepath,
            what=list(
              character(),
              integer(),
              NULL,
              integer()
              )
            )
  
  names(x)=c("chr","start","end","count")
  as.data.frame(x[c(1,2,4)],stringsAsFactors=F)
}
