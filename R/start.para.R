start.para <-
function(ncore, initExpr,...) {
  sfInit(parallel=TRUE, cpus=ncore, nostart=if.parallel(ncore),...)
  sfExport("initExpr")
  .dump <- sfClusterEval(eval(initExpr))
}
