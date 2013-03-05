start.para <-
function(ncore, varlist, type="PSOCK", ...) {
  if(ncore>1) cl <- makeCluster(getOption("cl.cores",ncore),
                                type=type,
                                ...)
  else cl <- NULL
  clusterExport(cl, varlist)
  .dump <- clusterEvalQ(cl, eval(dptws.initexpr))
  ##.dump <- sfClusterEval(eval(initExpr))
  cl
}
