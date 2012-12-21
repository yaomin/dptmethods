pattSummary.2dataframe <-
function(.out) {
  if(length(.out) ==1) {
    out <- data.frame(rbind(.out[[1]]$counts),
                      N=.out[[1]]$N,
                      best.cutoff=.out[[1]]$best.cutoff,
                      best.valid=.out[[1]]$best.valid,
                      best.pattern=.out[[1]]$best.pattern,
                      diff.score=.out[[1]]$diff.score)
                      
  } else {
    out <- data.frame(ldply(.out, function(x) x$counts)[,-1],
                    N=ldply(.out, function(x) x$N)[,-1],
                    best.cutoff=ldply(.out, function(x) x$best.cutoff)[,-1],
                    best.valid=ldply(.out, function(x) x$best.valid)[,-1],
                    best.pattern=ldply(.out, function(x) x$best.pattern)[,-1],
                    diff.score=ldply(.out, function(x) x$diff.score)[,-1])
  }
  names(out) <- c(names(out)[-(ncol(out)-(0:4))],
                  c("N","best.cutoff","best.valid","best.pattern","diff.score"))
  out
}
