node.set.significance.aux <-
function(anode, Pr.cutoff=0.05) {
  if(!is.leaf(anode)) {
    this.Pr <- attr(anode, 'Pr')
    ##if(attr(anode, 'Pr') < Pr.cutoff) attr(anode, 'Sig') <- TRUE else attr(anode, 'Sig') <- FALSE
    if(!is.na(this.Pr) && this.Pr < Pr.cutoff) attr(anode, 'Sig') <- TRUE
    else attr(anode, 'Sig') <- FALSE
  }
  else attr(anode, 'Sig') <- NA
  anode
}
