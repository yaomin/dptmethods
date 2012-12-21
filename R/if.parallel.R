if.parallel <-
function(ncpus) {
  if(length(ncpus)<1 || ncpus < 2) para.nostart=T else para.nostart=F
  para.nostart
}
