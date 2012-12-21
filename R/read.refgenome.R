read.refgenome <-
function(chrs, path="./FASTA/", ext=".fa") {
  cat("\n")
  reflist <- sapply(chrs, function(x) {cat(x);
                                       tst <- read.DNAStringSet(paste(path,x,ext,sep=""));
                                       cat("\tDONE\n");
                                       tst})
  names(reflist) <- chrs
  reflist
}
