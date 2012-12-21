map.regions <-
function(report.sig,
                        report,
                        chrom,
                        methylType,
                        patts=c("**1","*1*"),
                        cutoff=5,
                        winsize=100)
{  ## parallel the process on multiple processes
  trf.sig <- subset(report.sig,
                    subset= chr==chrom & MethylType==methylType)
  trf <- subset(report,
                subset= chr== chrom)

  reg.map <- map.regs.chr(report.sig,
                          report,
                          chrom,
                          methylType,
                          patts=patts,
                          cutoff=cutoff,
                          winsize=winsize)

  tst <- overlapProp.byID(trf.sig, trf, reg.map[,1:2], reg.map[,3:4])
  tst.subset <- overlap.commonSet(trf.sig, trf, reg.map[,1:2], reg.map[,3:4])

  tst1 <- overlapProp.byID(trf.sig, trf, reg.map[,1:2], reg.map[,5:6])
  tst1.subset <- overlap.commonSet(trf.sig, trf, reg.map[,1:2], reg.map[,5:6])

  tst.p2df <- data.frame(prop1=tst$p2, prop2=tst1$p2)

  tst.p2df.na <- apply(tst.p2df, 1, function(x) all(is.na(x)))

  tst.p2df.nona <- tst.p2df[!tst.p2df.na,]

  tst.p2df.choose <- rep(0, nrow(tst.p2df))

  tst.p2df.choose[!tst.p2df.na] <- apply(tst.p2df.nona, 1, which.max)


  tst.regions <- cbind(tst.subset$dat1,
                       (!tst.p2df.na),
                       tst.p2df,
                       tst.subset$dat2,
                       reg.map)

  tst.regions[tst.p2df.choose==0,] <- data.frame(
                tst.subset$dat1[tst.p2df.choose==0,],
                (!tst.p2df.na)[tst.p2df.choose==0],
                tst.p2df[tst.p2df.choose==0,],
                tst.subset$dat1[tst.p2df.choose==0,1:12],
                reg.map[tst.p2df.choose==0,]
                )
  tst.regions[tst.p2df.choose==1,] <- cbind(
                tst.subset$dat1[tst.p2df.choose==1,],
                (!tst.p2df.na)[tst.p2df.choose==1],
                tst.p2df[tst.p2df.choose==1,],
                tst.subset$dat2[tst.p2df.choose==1,],
                reg.map[tst.p2df.choose==1,])
  tst.regions[tst.p2df.choose==2,] <- cbind(
                tst.subset$dat1[tst.p2df.choose==2,],
                (!tst.p2df.na)[tst.p2df.choose==2],
                tst.p2df[tst.p2df.choose==2,],
                tst1.subset$dat2[tst.p2df.choose==2,],
                reg.map[tst.p2df.choose==2,])
  names(tst.regions) <- c(names(tst.subset$dat1),
                          "R.exist",
                          paste("R", names(tst.p2df), sep="."),
                          paste("R", names(tst.subset$dat2), sep="."),
                          names(reg.map))
  tst.regions
}
