start.para <-
function(ncore, initExpr, type="PSOCK", ...) {
  if(ncore>1) cl <- makeCluster(getOption("cl.cores",ncore),
                                type=type,
                                ...)
  else cl <- NULL
  .dump <- clusterEvalQ(cl, initExpr)
  ##.dump <- sfClusterEval(eval(initExpr))
  cl
}
