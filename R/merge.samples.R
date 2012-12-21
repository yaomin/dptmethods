#' @export
merge.samples <-
function(rootpath,mappedFiles,outpath) {
  nf <- length(mappedFiles)
  if(missing(outpath)) outpath <- get.ws.path("joinSample")

  n.reads <- vector("numeric", length=nf)
  count.tbs <- vector("list", length=nf)
  
  if(nf==1) {
    cat("loading aligned reads data...\n")
    in.a <- load.mappedFile(file.path(rootpath, mappedFiles))

    n.reads[1] <- get.n.reads.sample(in.a)
    count.tbs[[1]] <- get.count.tb.sample(in.a)
    
    chrs <- unique(in.a$chr)
    cat("merging samples...\n")
    out <- NULL
    for (chr in chrs) {
      cat("",tabsps(chr))
      chrom <- chr
      in.a.1 <- subset(in.a, subset=chr==chrom, select=c("start","count"))
      out <- c(out,list(R.merger(in.a.1)))
      cat("done\n")
    }
    names(out) <- chrs
  } else if(nf >1) {
    cat("loading first two aligned reads data...\n")
    in.a <- load.mappedFile(file.path(rootpath, mappedFiles[1]))
    in.b <- load.mappedFile(file.path(rootpath, mappedFiles[2]))

    n.reads[1] <- get.n.reads.sample(in.a)
    count.tbs[[1]] <- get.count.tb.sample(in.a)
    n.reads[2] <- get.n.reads.sample(in.b)
    count.tbs[[2]] <- get.count.tb.sample(in.b)
    
    chrs <- union(unique(in.a$chr), unique(in.b$chr))
    out <- NULL
    cat("merging first two samples...\n")
    for (chr in chrs) {
      cat("",tabsps(chr))
      chrom <- chr
      in.a.1 <- subset(in.a, subset=chr==chrom, select=c("start","count"))
      in.b.1 <- subset(in.b, subset=chr==chrom, select=c("start","count"))      
      out <- c(out, list(R.merger(in.a.1, in.b.1)))
      cat("DONE\n")
    }
    names(out) <- chrs
    if(nf > 2) {
      for(j in seq(3,nf)) {
        cat("loading aligned reads data",j,"...\n")
        in.j <- load.mappedFile(file.path(rootpath, mappedFiles[j]))

        n.reads[j] <- get.n.reads.sample(in.j)
        count.tbs[[j]] <- get.count.tb.sample(in.j)
        
        chrs <- union(chrs, unique(in.j$chr))
        cat("merging sample",j,"...\n")
        for (chr in chrs) {
          cat(" ", tabsps(chr))
          chrom <- chr
          in.j.1 <- subset(in.j, subset=chr==chrom, select=c("start","count"))
          if(!(chr %in% names(out)))
            out <- c(out, list(R.merger(subset(out[[names(out)[1]]], subset=NA), in.j.1)))
          else out[[chr]] <- R.merger(out[[chr]], in.j.1)
          cat("done\n")
        }
        names(out) <- chrs
      }
    }
  } else stop("mappedFiles not provided!")

  cat("writting merged data to workspace...\n")
  for (chr in names(out)) {
    ##browser()
    cat("",tabsps(chr))
    obj.name <- paste("s", chr, sep=".")
    assign(obj.name, as.data.frame(out[[chr]]))
    save.ws(ws="preprocess", what=obj.name,chr=chr, outpath=outpath)
  }
  cat("writting n.reads and count.table to workspace...\n")
  cat("",tabsps("n.reads"))
  save(n.reads, file=file.path(get.wsoption.path("joinSample"),"nreads.RData"))
  cat("done\n")
  cat("",tabsps("count.tbs"))
  save(count.tbs, file=file.path(get.wsoption.path("joinSample"),"countdistr.RData"))
  cat("done\n")
}
