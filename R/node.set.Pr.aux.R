node.set.Pr.aux <-
function(anode, n, p, k.cutoff=100,win.size=100) {
  if(!is.leaf(anode)) {
    value <- node.compute.Pr(anode, n, p, k.cutoff, win.size) }
  else value <- NA
  attr(anode, 'Pr') <- value
  anode
}
