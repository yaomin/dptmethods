#' @export
get.reads.fromJoinedSamples <-
function(inputPath, chr, threshold=1, pos=-1,
         pos.out=parent.frame(), aug.fun=max,
         sample.select=NULL) {
  dpt.options <- get("dpt.options", pos=pos)
  reads.datafile <- paste(inputPath,
                          paste(dpt.options[[c("output","joinSample")]],chr,"RData", sep="."),
                          sep="/")
  lf <- load(reads.datafile)

  s.cols <- grep("^s[0-9]*",names(get(lf)))
  if(!is.null(sample.select)) s.cols <- s.cols[sample.select]
  n.sample <- length(s.cols)

  # If length of threshold is 1, it will be applied to all samples
  # Otherwise, it is a vector of sample-specific cutoffs (window must pass at least 1 sample-specific cutoff to be retained)
  if(length(threshold)==1)
  {
    reads <- get(lf)[apply(get(lf)[,s.cols,drop=F], 1, aug.fun) > threshold,]
  } else {
    # Apply sample.select to threshold if it is a vector of per-sample cutoffs
    threshold <- threshold[sample.select]
    reads <- get(lf)[rowSums(get(lf)[,s.cols,drop=F] > threshold)>0,]
  }

  reads.poi <- reads$poi
  reads.names <- names(reads)[s.cols]
  assign("s.cols", s.cols, pos = pos.out)
  assign("reads", reads, pos=pos.out)
  assign("reads.poi", reads.poi, pos=pos.out)
  assign("reads.names", reads.names, pos=pos.out)
}
