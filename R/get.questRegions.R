get.questRegions <-
function(fasta.file) {
  qnames <- names(read.DNAStringSet(fasta.file))
  regions <- extract.info.fromQuESTRegioOutput(qnames)
  dlply(regions, 'chr', function(x) x)
}
