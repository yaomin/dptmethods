#'  @export

cleanup.on.gmeans <-
function(regions) {
  
  cat("clean up non-sharp sites ...\n")
  regions$patt <- regions$pattern
  uniquepatts <- as.character(unique(regions$patt))
  patt111 <- uniquepatts[sapply(strsplit(uniquepatts,
                                         split=""),
                                function(x) all(x=="1"))]
  ###... remaking the 0/1 decision based on LDA
  training <- subset(regions, subset=pattern!=patt111)
  lda.m <- with(training, lda(as.matrix(c(mean.0, mean.1)),
                              grouping=c(rep(0, length(mean.0)), rep(1, length(mean.1)))))
  ###... update site patterns based on LDA model
  regions.hyper.no111 <- subset(regions, subset=MethylType=="Hyper"&pattern!=patt111)
  sharp.hyper.no111 <- with(regions.hyper.no111,
                            predict(lda.m, as.matrix(mean.0))$class==0 & predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hyper.no111$pattern[!sharp.hyper.no111] <- patt111
  cat("num of validated sharp sites in hyper not-all-1 patterns :",
      sum(sharp.hyper.no111),
      "/",
      nrow(regions.hyper.no111) ,
      "\n")
  
  regions.hyper.111 <- subset(regions, subset=MethylType=="Hyper"&pattern==patt111)
  sharp.hyper.111 <- with(regions.hyper.111,
                          predict(lda.m, as.matrix(mean.0))$class==0 & predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hyper.111$pattern[sharp.hyper.111] <- regions.hyper.111$contrast[sharp.hyper.111]
  cat("num of validated sharp sites in hyper all-1 patterns :",
      sum(sharp.hyper.111),
      "/",
      nrow(regions.hyper.111) ,
      "\n")
   
  regions.hypo.no111 <- subset(regions, subset=MethylType=="Hypo"&pattern!=patt111)
  sharp.hypo.no111 <- with(regions.hypo.no111,
                           predict(lda.m, as.matrix(mean.0))$class==0&predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hypo.no111$pattern[!sharp.hypo.no111] <- patt111
  cat("num of validated sharp sites in hypo not-all-1 patterns :",
      sum(sharp.hypo.no111),
      "/",
      nrow(regions.hypo.no111) ,
      "\n")
  

  regions.hypo.111 <- subset(regions, subset=MethylType=="Hypo"&pattern==patt111)
  sharp.hypo.111 <- with(regions.hypo.111,
                         predict(lda.m, as.matrix(mean.0))$class==0&predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hypo.111.patts <- regions.hypo.111$contrast[sharp.hypo.111]
  regions.hypo.111.patterns <- gsub("b","0", gsub("0", "1", gsub("1","b", regions.hypo.111.patts)))
  regions.hypo.111$pattern[sharp.hypo.111] <- regions.hypo.111.patterns
  cat("num of validated sharp sites in hypo all-1 patterns :",
      sum(sharp.hypo.111),
      "/",
      nrow(regions.hypo.111) ,
      "\n")
  
  rbind(regions.hyper.no111, regions.hyper.111, regions.hypo.no111, regions.hypo.111)
}
