cutoff.sites <-
function(q.dt,
                         qvalue.cutoff=NA, value.cutoff=NA,
                         value.var="value.1", qvalue.var="qvalue") {
  ## remove rows with NAs in either value or qvalue
  to.rm <- c(which(is.na(q.dt[[qvalue.var]])),
             which(is.na(q.dt[[value.var]])))
  if(length(to.rm) > 0) {
    q.dt <- q.dt[-to.rm,]
  }
  
  if(is.na(qvalue.cutoff) & is.na(value.cutoff)) {
  } else if(is.na(qvalue.cutoff)) {
    q.dt <- q.dt[abs(q.dt[[value.var]]) > value.cutoff,]
    ## value
  } else if(is.na(value.cutoff)) {
    q.dt <- q.dt[q.dt[[qvalue.var]]<qvalue.cutoff,]
  } else {
    q.dt <- q.dt[abs(q.dt[[value.var]]) > value.cutoff & q.dt[[qvalue.var]]<qvalue.cutoff,]
  }
  q.dt
}
