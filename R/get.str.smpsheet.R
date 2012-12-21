get.str.smpsheet <-
function(smpsheetStr) {
  smpsheet <- gsub('-','',grep('samplesheet', commandArgs(), v=T))
  if(missing(smpsheetStr)) smpsheet else smpsheetStr
}
