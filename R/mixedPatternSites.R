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
      idx.mixed <- sites[[i]]%over%sites[[j]]
      mix.j <- sites[[i]][idx.mixed]
      mix.i <- if(is.null(mix.i)) mix.j else c(mix.i, mix.j)
    }
    mix.i <- unique(mix.i)
    mixed.r <- if(is.null(mixed.r)) mix.i else c(mixed.r, mix.i)
  }
  unique(mixed.r)
}

reAssign.mixedPatternSites <- function(mixedSites, e.TF) {

  mixedSites.disjoined <- disjoin(mixedSites)
  e.TF.r <- ranges(e.TF)[[1]]
  ## remove P000 on two ends
  e.sub.0 <- e.TF[e.TF.r %over% mixedSites.disjoined,]
  e.sub.v <- as.data.frame(values(e.sub.0))
  idx000 <- (apply(e.sub.v[,-(1:2)], 1, sum)==0)&e.sub.v[grep("P0+$", names(e.sub.v))]
  sub000 <- ranges(e.sub.0)[[1]][idx000]
  site000 <- reduce(sub000)

  mtch1 <- subsetByOverlaps(site000, mixedSites.disjoined, type="start")
  mtch2 <- subsetByOverlaps(site000, mixedSites.disjoined, type="end")
  mtched <- c(mtch1, mtch2)

  .disjoin <- disjoin(c(mixedSites.disjoined, mtched))
  .disjoin2 <- .disjoin[!(.disjoin%over%mtched)]

  ## compuate the best pattern for each disjoined site
  mtch.idx2 <- findOverlaps(e.TF.r, .disjoin2)

  .summary <-  by(as.data.frame(values(e.TF[queryHits(mtch.idx2),]))[,-1],
                  list(subjectHits(mtch.idx2)),
                  function(x, cutoff) patt.summary(x, cutoff),
                  cutoff=0.5,
                  simplify=F)

  .summary.df <- pattSummary.2dataframe(.summary)

  out.mixed <- split(.disjoin2[unique(subjectHits(mtch.idx2)),], .summary.df$best.pattern)
  out.mixed[-grep("P0+$", names(out.mixed))]
}

handle.mixedPatternSites <- function(sites, e.TF,
                                     mixed.filter.cut =120,
                                     join.gap=NULL) {
  if(is.null(join.gap)) join.gap <- mixed.filter.cut
  if(is(sites, "list")) sites <- RangesList(sites)
  mixed.r <- find.mixedPatternSites(sites, mixed.filter.cut)
  if(length(mixed.r)!=0)
  {
  	mixed.sites <- reAssign.mixedPatternSites(mixed.r, e.TF)
  	##pure.sites <- sites[!(sites%in%mixed.sites)]
  	pure.sites <- IRangesList(lapply(names(sites),
                                   function(x) sites[[x]][!(sites[[x]]%over%mixed.r)]))
  	names(pure.sites) <- names(sites)
  	IRangesList(lapply(combine.2rl(pure.sites, mixed.sites),
                     reduce,
                     min.gapwidth=join.gap))
  } else
  {
  	sites
  }
}
