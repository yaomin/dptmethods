find.cutoff <-
function(cutoff, paras, range=c(-10,10)) {
  ret <- optimize(find.cutoff.objective,
                  range,
                  paras=paras,
                  cutoff=cutoff)
  round(ret$minimum, 3)
}
