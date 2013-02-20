batch.jobs <-
function(caller, nc, nb, chrs=NULL, by.chr=T, retry.n=3, waiting=90) {
  browser()
  .failed <- NULL
  ncore <- nc
  if(!by.chr) {
    .failed <- eval(caller)
  } else {
    .continue <- T
    .ntry <- 1
    while(.continue) {
      if(.ntry>1) cat("try = ", .ntry, "\n")
      .jobs <- (seq(along=chrs)-1)%/%nb
      .failed.try <- NULL
      for (i in unique(.jobs)) {
        .chrs <- chrs[.jobs==i]
        for (chr in .chrs){
          multicore::mcparallel(caller)
          dpt.syshold(waiting)
        }
        .fl <- c(unlist(multicore::collect()))
        .failed.try <- c(.failed.try, .fl)
      }
      chrs <- .failed.try
      .ntry <- .ntry+1
      .continue <- (.ntry<=retry.n)&(length(.failed.try)>0)
    }
    .failed <- .failed.try
  }
  return(.failed)
}
