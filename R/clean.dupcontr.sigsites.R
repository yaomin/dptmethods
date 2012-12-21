clean.dupcontr.sigsites <-
function(sigsites,
                                    p.value.col)
{
  ## this is a inefficient implementation, optimize it whenever possible
  out.1 <- by(sigsites,
              list(chr.ID=paste(sigsites$chr,
                     sigsites$pattern,
                     sigsites$ID,                     
                     sep=":")),
              function(x) {
                nr <- dim(x)[1]
                if(nr<2) x
                else x[which.min(x[,p.value.col]),]
              })
  out.2 <- ldply(out.1, function(x) x)
  out.2
}
