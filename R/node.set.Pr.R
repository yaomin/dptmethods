node.set.Pr <-
function(anode,n,p,k.cutoff=100, win.size=100) {
  dendrapply(anode, node.set.Pr.aux, n=n, p=p, k.cutoff=k.cutoff, win.size=win.size)
}
