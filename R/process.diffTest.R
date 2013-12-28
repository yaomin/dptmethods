#' @export
process.diffTest <-
function(cl,sites.js.ex, events.wins, events.ctr, testmode.rep,
                             pmeans.norm, pmeans.var,
                             ncore,
                             sample.label,
                             subject.label) {
  patts <- unique(sites.js.ex$pattern)
  events.wins.matched <- events.wins[,match.patt.pattmx(patts, events.wins), drop=F]
  events.ctr.matched <- events.ctr[match.patt.pattmx(patts, events.wins)]
  patt.n <- length(patts)
  wins.test <- vector("list", length=patt.n)
  cat("testing ...\n")
  if(testmode.rep) {
    for (i in seq(patt.n)) {
      wins.test[[i]] <- list()
      for(j in seq(length(events.ctr.matched[[i]]))) {
        cat(as.character(patts[i]), "\t", j, "\n")
        wins.test[[i]][[j]] <- test.effect.rep(cl,
                                               subset(sites.js.ex, subset=pattern==patts[i]),
                                               pmeans.norm,
                                               pmeans.var,
                                               contr= events.ctr.matched[[i]][[j]],
                                               n.core = ncore,
                                               group.label=sample.label,
                                               rep.label=subject.label)
        attr(wins.test[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
        attr(wins.test[[i]][[j]],'pattern') <- events.wins.matched[,i]
      }
    }
  } else {
    for (i in seq(patt.n)) {
      wins.test[[i]] <- list()
      for(j in seq(length(events.ctr.matched[[i]]))) {
        cat(as.character(patts[i]), "\t", j, "\n")
        wins.test[[i]][[j]] <- test.effect(cl,
                                           subset(sites.js.ex, subset=pattern==patts[i]),
                                           pmeans.norm,
                                           pmeans.var,
                                           contr= events.ctr.matched[[i]][[j]],
                                           n.core = ncore, group.label=sample.label)
        attr(wins.test[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
        attr(wins.test[[i]][[j]],'pattern') <- events.wins.matched[,i]
      }
    }
  }

  cat("report:\n")
  wins.test.report <- vector("list", length=patt.n)
  for (i in seq(patt.n)) {
    wins.test.report[[i]] <- list()
    for(j in seq(length(events.ctr.matched[[i]]))) {
      cat(as.character(patts[i]), "\t", j, "\n")
      wins.test.report[[i]][[j]] <- combine.report(subset(sites.js.ex, subset=pattern==patts[i]),
                                                   wins.test[[i]][[j]])
      attr(wins.test.report[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
      attr(wins.test.report[[i]][[j]],'pattern') <- events.wins.matched[,i]
    }
  }
  list(wins.test=wins.test, wins.test.report=wins.test.report)
}
