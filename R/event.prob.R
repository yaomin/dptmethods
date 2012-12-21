event.prob <-
function(x, sample.label,event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), k=1e-36,
                       which.score=.25, npart=10,
                       para.mode=1)
  ## assume x is in w3.res format as above
  ## p is the p3
{
  ##browser()
  escore <- event.prob.mx(x, event, p, k)
  
  if(para.mode==1) {
    output <- paraApply(escore, 1,
                        function(x, which.score,smp.label) escore.comp(x, sample.label,which.score),
                        which.score=which.score,
                        ncore=npart,
                        mode=real64(),
                        mode.size=8,
                        smp.label=sample.label)
  } else {
    output <- paraApply2(escore,
                         1,
                         function(x, which.score,smp.label) escore.comp(x, sample.label,which.score),
                         which.score=which.score,
                         ncore=npart,
                         smp.label=sample.label)
  }
  output
}
