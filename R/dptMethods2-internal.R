.dump <-
TRUE
.First <-
function(){
  vflag <- F
  require(MCMCpack,quietly=vflag)
  require(snowfall,quietly=vflag)
  require(nlme,quietly=vflag)
  require(gmodels,quietly=vflag)
  require(MASS,quietly=vflag)
  require(plyr,quietly=vflag)
  require(preprocessCore,quietly=vflag)
  require(inline,quietly=vflag)
  require(Rcpp,quietly=vflag)
  require(IRanges,quietly=vflag)
  require(Biostrings,quietly=vflag)
  require(mmap,quietly=vflag)
  require(qvalue, quietly=vflag)
}
