#'  @export
compute.pmeans <-
function(res, reads.poi, reads.names, count.table, pos.out=parent.frame()) {

  ## if(!missing(ws.mixPois.path)) load.mixPoiData(ws.mixPois.path)
  ## res <- res[!unlist(lapply(res, is.null))]
  
  pmean.res <- as.data.frame(lapply(res, get.postMean, log=F))
  row.names(pmean.res) <- reads.poi
  names(pmean.res) <- reads.names

  pmean.var <- as.data.frame(lapply(res, get.postmean.var))
  row.names(pmean.var) <- reads.poi
  names(pmean.var) <- reads.names

  pmeans <- list(##lambda2 = lambda2.res,
                 mean = pmean.res,
                 var = pmean.var)
  
  ##pmeans.norm <- normalize.sig(pmeans$mean, pmeans$lambda2)
  ##pmeans.norm <- normalize.pmeans(pmeans$mean, countab)
  pmeans.norm <- as.data.frame(normalize.quantiles(as.matrix(pmean.res)))
  row.names(pmeans.norm) <- reads.poi
  names(pmeans.norm) <- reads.names

  ### bg corrected
  pmeans.bgcorrect.norm <- as.data.frame(normalize.quantiles(as.matrix(normalize.pmeans(pmean.res,
                                                                                        count.table,
                                                                                        bg.correct=T,
                                                                                        out.int=F))))
  row.names(pmeans.bgcorrect.norm) <- reads.poi
  names(pmeans.bgcorrect.norm) <- reads.names

  ### bg normalized but not corrected
  pmeans.bgnorm.norm <- as.data.frame(normalize.quantiles(as.matrix(normalize.pmeans(pmean.res,
                                                                                        count.table,
                                                                                        bg.correct=F,
                                                                                        out.int=F))))
  row.names(pmeans.bgnorm.norm) <- reads.poi
  names(pmeans.bgnorm.norm) <- reads.names

  assign("pmeans", pmeans, pos = pos.out)
  assign("pmeans.norm", pmeans.norm, pos = pos.out)
  assign("pmeans.bgnorm.norm", pmeans.bgnorm.norm, pos = pos.out)
  assign("pmeans.bgcorrect.norm", pmeans.bgcorrect.norm, pos = pos.out)
  
}
