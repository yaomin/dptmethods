get.n.cpus <-
function(n.cpus) {
  ncpus <- as.numeric(gsub("cpus","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) ==0) ncpus <- as.numeric(gsub("ncpus","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) ==0) ncpus <- as.numeric(gsub("ncores","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) == 0) ncpus <- 1
  if(missing(n.cpus)) ncpus else n.cpus
}
