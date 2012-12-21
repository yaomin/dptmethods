comb.hyperHypo <-
function(wins.test.hyper,
                           wins.test.hypo,
                           v.varname,
                           p.varname) {
  wins.test.comb <- rbind(subset(wins.test.hyper, subset=wins.test.hyper[, v.varname] > 0),
                          subset(wins.test.hypo, subset=wins.test.hypo[,v.varname] < 0))
  wins.test.cutoff <- select.wins.test.element(wins.test.comb,
                                               pvalue.colname=p.varname,
                                               effect.colname=v.varname)
  wins.test.cutoff
}
