## Filename: mixedPatternSites.R
## By: Yaomin Xu (yaomin@gmail.com)
## Date: Tue Mar 26 20:44:19 2013
## Notes:   handle the mixed pattern sites
## =============================================================================

find.mixedPatternSites <- function(sites, width.cutoff=120) {
  if(is(sites, "list")) sites <- RangesList(sites)
  sites <- sites[width(sites)>width.cutoff]
  patts <- names(sites)
  idx.n <- seq(length(sites))
  ## loop
  mixed.r <- NULL
  for (i in idx.n) {
    mix.i <- NULL
    cat(patts[i], ":\n")
    for (j in setdiff(idx.n, i)) {
      cat("\tcompare to", patts[j],"\n")
      idx.mixed <- sites[[i]]%in%sites[[j]]
      mix.j <- sites[[i]][idx.mixed]
      mix.i <- if(is.null(mix.i)) mix.j else c(mix.i, mix.j)
    }
    mix.i <- unique(mix.i)
    mixed.r <- if(is.null(mixed.r)) mix.i else c(mixed.r, mix.i)
  }
  unique(mixed.r)
}

reAssign.mixedPatternSites <- function(mixedSites, e.TF, cutoff=0.5) {
  
  mixedSites.disjoined <- disjoin(mixedSites)
  e.sub.0 <- e.TF[ranges(e.TF)[[1]] %in% mixedSites.disjoined,]
  e.sub.v <- as.data.frame(values(e.sub.0))
  idx000 <- (apply(e.sub.v[,-(1:2)], 1, sum)==0)&e.sub.v$P000
  sub000 <- ranges(e.sub.0)[[1]][idx000]
  site000 <- reduce(sub000)

  mtch1 <- subsetByOverlaps(site000, mixedSites.disjoined, type="start")
  mtch2 <- subsetByOverlaps(site000, mixedSites.disjoined, type="end")
  mtched <- c(mtch1, mtch2)

  .disjoin <- disjoin(c(mixedSites.disjoined, mtched))
  .disjoin2 <- .disjoin[!(.disjoin%in%mtched)]

  e.sub <- e.TF[ranges(e.TF)[[1]] %in% .disjoin2,]
  e.ranges <- ranges(e.sub)[[1]]

  .disjoin3 <- .disjoin2[.disjoin2%in%e.ranges]

  e.sub2 <- e.TF[ranges(e.TF)[[1]] %in% .disjoin3,]
  e.ranges2 <- ranges(e.sub2)[[1]]
  
  .fac <- rep(NA, length(e.sub2[[1]]))
  
  .fac <- match(e.ranges2, .disjoin3)
  
  .summary <-  by(as.data.frame(values(e.sub2))[,-1],
                  list(.fac),
                  function(x, cutoff) patt.summary(x, cutoff),
                  cutoff=cutoff,
                  simplify=F)
  .summary.df <- pattSummary.2dataframe(.summary)

  out.mixed <- split(.disjoin3, .summary.df$best.pattern)[-1] ## remove P000
  out.mixed
}

handle.mixedPatternSites <- function(sites, e.TF,
                                     mixed.filter.cut =120,
                                     cutoff=0.5) {
  if(is(sites, "list")) sites <- RangesList(sites)
  mixed.r <- find.mixedPatternSites(sites, mixed.filter.cut)
  mixed.sites <- reAssign.mixedPatternSites(mixed.r, e.TF)
  ##pure.sites <- sites[!(sites%in%mixed.sites)]
  pure.sites <- IRangesList(lapply(names(sites),
                                   function(x) sites[[x]][!(sites[[x]]%in%mixed.r)]))
  names(pure.sites) <- names(sites)
  combine.2rl(pure.sites, mixed.sites)
}
