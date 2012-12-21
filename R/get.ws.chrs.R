#' @export
get.ws.chrs <-
function(which=c("mixPoisson",
                          "pattRecog",
                          "diffTest",
                          "joinSample",
                          "full"
                          ),
                        pos=-1
                        )
{
  which <- match.arg(which)
  dpt.options <- get("dpt.options", pos)
  if(which != "full") {
    files <- list.files(get.ws.path(which, pos=pos))
    if(which == "mixPoisson") {
      filestrs <- grep(paste(dpt.options[[c('output','mixPoisson')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "pattRecog") {
      filestrs <- grep(paste(dpt.options[[c('output','pattRecog')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "diffTest") {
      filestrs <- grep(paste(dpt.options[[c('output','diffTest')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "joinSample") {
      filestrs <- grep(paste(dpt.options[[c('output','joinSample')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which =="full") {
      filestrs <- NULL
    } else stop("Not a valid option: ", which)
    out <- sapply(filestrs,
                  function(x) {
                    this.str <- unlist((strsplit(x, split="\\.")))
                    this.str[length(this.str)-1]
                  }
                  )
  }
  if(which == "full") out <- c(
       "chr10","chr11","chr12","chr13","chr14",
       "chr15","chr16","chr17","chr18","chr19",
       "chr1","chr20","chr21","chr22","chr2",
       "chr3","chr4","chr5","chr6","chr7","chr8",
       "chr9","chrX","chrY"
       )
  as.vector(out)
}
