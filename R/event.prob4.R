event.prob4 <-
function(x, sample.label,cl,
         event=c(0,0,0,1,1,1,1,1,1), e.cut=0.5,p=rep(0.33,9),
                       which.score=.25, npart=10,
                       para.mode=1,p.fun="median")
  ## assume x is in w3.res format as above
  ## p is the p3
{

  escore <- event.prob.mx3(x, event, p)
  
  if(para.mode==1) {
      .P <- paraApply(cl,
		      escore$x, 1,
		      function(x, which.score,smp.label) escore.comp3(x, sample.label,which.score),
		      which.score=which.score,
		      ncore=npart,
		      mode=real64(),
		      mode.size=8,
		      smp.label=sample.label)
  } else {
      .P <- paraApply2(cl,
		       escore$x,
		       1,
		       function(y, which.score,smp.label) escore.comp3(y, sample.label,which.score),
		       which.score=which.score,
		       ncore=npart,
		       smp.label=sample.label)
  }
  events.idx <- apply(.P > e.cut, 1, all)
  .P[!events.idx,] <- 0
  .p <- tapply(p, sample.label, p.fun)

  out <- apply(t(apply(.P,1, function(x) x*.p)), 1, prod)

  out
}
