escore.comp3 <-
function(x, group, which.p) {
    grp <- as.factor(group)
    p.len <- length(which.p)
    if(p.len>1) {
	stopifnot(nlevels(grp)!=p.len)
	sapply(seq(p.len),
	       function(i) apply(x[,grp==levels(grp)[i]],
				 1,
				 quantile,
				 probs=which.p[i]))
    } else if(p.len==1) {
	t(apply(x, 1, function(y) {
	      tapply(y,
		     group,
		     function(z, which.p) quantile(z, p=which.p),
		     which.p=which.p)}))
    } else stop("wrong length for which.p")
}
