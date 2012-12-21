node.set.significance <-
function(anode, Pr.cutoff=0.001) {
  dendrapply(anode, node.set.significance.aux,Pr.cutoff=Pr.cutoff)
}
