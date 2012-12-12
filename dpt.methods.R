## Filename: DiffSeq.methods.R
## By: Yaomin Xu
## Date: Sep 11, 2010
## Notes: DiffSeq R methods
## =============================================================================

###. merge samples: newer with faster implementation and less memory demand
##___________________________________________________________________
load.mappedFile <- function(filepath, cat.pre="", cat.after="\t") {
  cat(cat.pre, basename(filepath),cat.after)
  x <- scan(filepath,
            what=list(
              character(),
              integer(),
              NULL,
              integer()
              )
            )
  
  names(x)=c("chr","start","end","count")
  as.data.frame(x[c(1,2,4)],stringsAsFactors=F)
}

tabsps <- function(str) {
  paste(str, paste(rep("\t", (nchar(str)%/%6<1)+1), collapse="", sep=""), sep="")
}

get.n.reads.sample <- function(mappedData) {
  sum(mappedData[,"count"])
}

get.count.tb.sample <- function(mappedData) {
  tb.df <- as.data.frame(table(mappedData[,"count"]))
  names(tb.df) <- c("count", "freq")
  tb.df
}

merge.samples <- function(rootpath,mappedFiles,outpath) {
  nf <- length(mappedFiles)
  if(missing(outpath)) outpath <- get.ws.path("joinSample")

  n.reads <- vector("numeric", length=nf)
  count.tbs <- vector("list", length=nf)
  
  if(nf==1) {
    cat("loading aligned reads data...\n")
    in.a <- load.mappedFile(file.path(rootpath, mappedFiles))

    n.reads[1] <- get.n.reads.sample(in.a)
    count.tbs[[1]] <- get.count.tb.sample(in.a)
    
    chrs <- unique(in.a$chr)
    cat("merging samples...\n")
    out <- NULL
    for (chr in chrs) {
      cat("",tabsps(chr))
      chrom <- chr
      in.a.1 <- subset(in.a, subset=chr==chrom, select=c("start","count"))
      out <- c(out,list(R.merger(in.a.1)))
      cat("done\n")
    }
    names(out) <- chrs
  } else if(nf >1) {
    cat("loading first two aligned reads data...\n")
    in.a <- load.mappedFile(file.path(rootpath, mappedFiles[1]))
    in.b <- load.mappedFile(file.path(rootpath, mappedFiles[2]))

    n.reads[1] <- get.n.reads.sample(in.a)
    count.tbs[[1]] <- get.count.tb.sample(in.a)
    n.reads[2] <- get.n.reads.sample(in.b)
    count.tbs[[2]] <- get.count.tb.sample(in.b)
    
    chrs <- union(unique(in.a$chr), unique(in.b$chr))
    out <- NULL
    cat("merging first two samples...\n")
    for (chr in chrs) {
      cat("",tabsps(chr))
      chrom <- chr
      in.a.1 <- subset(in.a, subset=chr==chrom, select=c("start","count"))
      in.b.1 <- subset(in.b, subset=chr==chrom, select=c("start","count"))      
      out <- c(out, list(R.merger(in.a.1, in.b.1)))
      cat("DONE\n")
    }
    names(out) <- chrs
    if(nf > 2) {
      for(j in seq(3,nf)) {
        cat("loading aligned reads data",j,"...\n")
        in.j <- load.mappedFile(file.path(rootpath, mappedFiles[j]))

        n.reads[j] <- get.n.reads.sample(in.j)
        count.tbs[[j]] <- get.count.tb.sample(in.j)
        
        chrs <- union(chrs, unique(in.j$chr))
        cat("merging sample",j,"...\n")
        for (chr in chrs) {
          cat(" ", tabsps(chr))
          chrom <- chr
          in.j.1 <- subset(in.j, subset=chr==chrom, select=c("start","count"))
          if(!(chr %in% names(out)))
            out <- c(out, list(R.merger(subset(out[[names(out)[1]]], subset=NA), in.j.1)))
          else out[[chr]] <- R.merger(out[[chr]], in.j.1)
          cat("done\n")
        }
        names(out) <- chrs
      }
    }
  } else stop("mappedFiles not provided!")

  cat("writting merged data to workspace...\n")
  for (chr in names(out)) {
    ##browser()
    cat("",tabsps(chr))
    obj.name <- paste("s", chr, sep=".")
    assign(obj.name, as.data.frame(out[[chr]]))
    save.ws(ws="preprocess", what=obj.name,chr=chr, outpath=outpath)
  }
  cat("writting n.reads and count.table to workspace...\n")
  cat("",tabsps("n.reads"))
  save(n.reads, file=file.path(get.wsoption.path("joinSample"),"nreads.RData"))
  cat("done\n")
  cat("",tabsps("count.tbs"))
  save(count.tbs, file=file.path(get.wsoption.path("joinSample"),"countdistr.RData"))
  cat("done\n")
}

R.merger <- function(smp1,smp2) {
  if(missing(smp2)) {
    out <- smp1
    dimnames(out)[[2]] <- c("poi", "s1")
  } else {
    poi <- union(smp1[,1], smp2[,1])
    x1 <- matrix(0, nrow=length(poi), ncol=ncol(smp1)-1)
    x2 <- rep(0, length(poi))
    x1[match(smp1[,1], poi),] <- smp1[,-1]
    x2[match(smp2[,1], poi)] <- smp2[,2]
    out <- cbind(poi, x1, x2)
    dimnames(out)[[2]] <- c("poi", paste("s", 1:(ncol(x1)+1), sep=""))
  }
  out
}

Rcpp.merger <- function(ids, smp) {
  src <-'
     string s1, s2;
     Rcpp::CharacterVector xid(id);
     Rcpp::NumericVector xpoi(poi);
     Rcpp::CharacterVector xsid(sid);
     Rcpp::NumericVector xspoi(spoi);
     Rcpp::NumericVector xscnt(scnt);
     int n_xid = xid.size();
     int n_xsid = xsid.size();
     Rcpp::NumericVector xcnt(n_xid);
     for (int i=0; i<n_xid; i++)
        for (int j=0; j<n_xsid; j++) {
          s1 = xid[i];
          s2 = xsid[j];
          if (s1 == s2 & xpoi[i]- xspoi[j] <1) xcnt[i] = xscnt[j];
          else xcnt[i] =0;
        };  
     return xcnt;
  '

  .merger <- cxxfunction(signature(id="character",
                                   poi="numeric",
                                   sid="character",
                                   spoi="numeric",
                                   scnt="numeric"),
                         includes=c("#include <string>", "using namespace std;"),
                         src,
                         plugin="Rcpp")

  .merger(ids$chr, ids$start, smp$chr, smp$start, smp$count)
}



###. join samples
##____________________________________________________________________

win.movingCount <- function(ares, win.size=4) {
  pulse.n <- dim(ares)[1]
  win.n <- pulse.n - win.size + 1
  res <- sapply(seq(win.n),
                function(x, dat, win.size) sum(dat[x - 1 + seq(win.size),2]),
                dat = ares,
                win.size=win.size)
  cbind(V2=ares[seq(win.n),1], V4=res)
}

truncate4wd <- function(x, n=13) {
  x[seq(2^n),]
}

###.. Extend samples
##--------------------------------------
extend.seq <- function(seqdat, start=NULL, end=NULL, by=100) {
  
  if(is.null(start)) if(nrow(seqdat) > 0) NA else start=min(seqdat[,1])
  if(is.null(end)) if(nrow(seqdat) > 0) NA else end = max(seqdat[,1])
  V2 = seq(from=start, to=end, by = by)
  V4 = rep(0, length(V2))
  exted <- as.data.frame(cbind(V2,
                               V4))

  if(nrow(seqdat) > 0) exted[match(seqdat[,1],exted$V2),2] <- seqdat[,2]
  exted
}


###.. Join samples (2 or more)
##--------------------------------------

range.poi <- function(s1, s2, which.col=2) {
  if(!missing(s2)) {
    if(nrow(s1) <1 & nrow(s2) <1) return(NA)
    else {
      range.s1 <- if(nrow(s1)==0) NA else range(s1[,which.col])
      range.s2 <- if(nrow(s2)==0) NA else range(s2[,which.col])
      range(range.s1, range.s2, na.rm=T)
    }
  }
  else {
    if(nrow(s1) <1) return(NA)
    else {
      range.s1 <- if(nrow(s1)==0) NA else range(s1[,which.col])
      range.s1
    }
  }
}

range.pois <- function(x, ..., which.col=2) {

  if(nargs()<4) range.poi(x,..., which.col=which.col)
  else {
    range(if(nrow(x) >0 ) range(x[,which.col]) else NA,
          range.pois(..., which.col=which.col), na.rm=T)
  }
}

join.samples <- function(x,...,by=100) {
  
  ends <- range.pois(x, ..., which.col=2)
  x.ext <- extend.seq(x[,c(2,4)], start=ends[1], end=ends[2],by=by)
  if(nargs()>2) {
    res <- lapply(list(...),function(x, ends) extend.seq(x[,c(2,4)],
                                                         start=ends[1],
                                                         end=ends[2],
                                                         by=by),
                  ends=ends)
    res.cnt <- as.data.frame(sapply(res, function(x) x[,2]))
    res.cmb <- cbind(x.ext, res.cnt)
  } else {
    res.cmb <- x.ext
  }
  names(res.cmb) <- c('poi', paste('s', seq(ncol(res.cmb)-1), sep=''))
  res.cmb
  
}

join.counts <- function(x, ...,by=100) {
  
  ends <- range.pois(x, ..., which.col=1)
  x.ext <- extend.seq(x, start=ends[1], end = ends[2],by=by)
  res <- lapply(list(...), function(x, ends) extend.seq(x,
                                                        start=ends[1],
                                                        end=ends[2],
                                                        by=by),
                ends = ends)
  res.cnt <- as.data.frame(sapply(res, function(x) x[,2]))
  res.cmb <- cbind(x.ext, res.cnt)
  names(res.cmb) <- c('poi', paste('c', seq(ncol(res.cmb)-1), sep=''))
  res.cmb
  
}

add.counts <- function(x, ...) {
  ## if(nargs() < 3) add.counts.2smp.v2(x, ...)
  ## else add.counts.2smp.v2(x, add.counts.2smp.v2(...))
  res <- join.samples(x, ...)
  res <- as.data.frame(cbind(res[,1], rowSums(res[,-1])))
  names(res) <- c('poi', 'c1')
  res
}

join.samples.chr <- function(chr, x, ..., by=100,filePrefix=NULL, pos=-1) {
  dpt.options <- get("dpt.options", pos=pos)
  if(is.null(filePrefix)) filePrefix <- paste(get.wsoption.path("joinSample", pos=pos),
                                              dpt.options[[c("output","joinSample")]],
                                              sep="/")
  arglist <- lapply(list(x,...), function(x) subset(x, subset=(V1 == chr)))
  filename <- paste(filePrefix, chr, "RData",sep=".")
  assign(paste("s",chr, sep="."), do.call(join.samples, args=c(arglist,by=by)), pos=pos)
  ##save.env <- if(pos==-1) globalenv() else as.environment(pos)
  save(list=(paste("s",chr, sep=".")), file=filename)
}

join.counts.chr <- function(chr, x, ...,by=100) {
  
  ## notice: the data structures of x,... are differet from those in join.samples.chr
  arglist <- lapply(list(x,...), function(x) x[[chr]])
  do.call(join.counts, args=c(arglist,by=by))
}

add.counts.chr <- function(chr, x, ...,by=100) {
  arglist <- lapply(list(x,...), function(x) subset(x, subset=(V1 == chr)))
  do.call(add.counts, args=c(arglist, by=by))
}
  

para.joinSamples <- function(x, ..., n.clusters=2, chrs=NULL, by=100, filePrefix=NULL, pos=-1) {
  dpt.options <- get("dpt.options", pos=pos)
  if(is.null(filePrefix)) filePrefix <- paste(get.wsoption.path("joinSample", pos=pos),
                                              dpt.options[[c("output","joinSample")]],
                                              sep="/")
  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  function(x) as.character(unique(x$V1)))))
  }
  
  res <- sfLapply(chrs, join.samples.chr, x, ..., filePrefix=filePrefix,by=by)
}

##TODO: Update snowfall handler. At least put it to ther outer main session
para.joinCounts <- function(x, ..., n.clusters=2, chrs=NULL,
                            sourcefile=NULL,
                            winsize=100) {
  
  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  names)))
  }
  
  sfInit(parallel = T, cpus = n.clusters)
  sfExportAll()
  ##sfSource(sourcefile)
  res <- sfLapply(chrs, join.counts.chr, x, ...,by=winsize)
  sfStop()
  
  names(res) <- chrs
  res
}



para.addCounts <- function(x, ..., n.clusters=2, chrs=NULL,
                           sourcefile=NULL) {

  if(is.null(chrs)) {
    chrs <- unique(unlist(lapply(list(x,...),  function(x) as.character(unique(x$V1)))))
    ##chrs <- unique(unlist(lapply(list(x,...),  names)))
  }
  
  sfInit(parallel = T, cpus = n.clusters)
  sfExportAll()
  ##sfSource(sourcefile)
  
  res <- sfLapply(chrs, add.counts.chr, x, ...)
  sfStop()
  names(res) <- chrs
  res
}

test.ratios <- function(dat, x.n, y.n, w.n=2500000){
 ## Assume we have data from join.samples
 ##require(rateratio.test)

 res <- apply(dat,
              1,
              function(adat, x.n,y.n, w.n)
              as.numeric(unlist(rateratio.test(adat, c(x.n, y.n)/w.n))[c(1,2,3,4,6,7)]),
              x.n=x.n,
              y.n=y.n,
              w.n=w.n)
 t(res)
}


test.ratios.chr <- function(chr, dat, x.n, y.n, w.n=2500000, cols=2:3) {
  ## test only on one chr
  
  dat.sub <- dat[[chr]][,cols]
  res <- as.data.frame(test.ratios(dat.sub, x.n, y.n, w.n))
  names(res) <- c('pval','ratio','rate1', 'rate2', 'ciL','ciU')
  res
}


para.testRatios <- function(dat, x.n, y.n, w.n=2500000, cols=2:3,
                            chrs=NULL, n.clusters=2,
                            sourcefile=NULL) {
  ## dat is assumed from para.joinSamples
  dat <- lapply(dat, consolidate.lows)
  if(is.null(chrs)) chrs <- seq(length(dat))
  else if(is.character(chrs)) chrs <- match(chrs, names(dat))
  sfInit(parallel = T, cpus = n.clusters)
  sfExport('dat', local=T)
  ##sfSource(sourcefile)
  res <- sfLapply(chrs, test.ratios.chr, dat=dat, x.n=x.n, y.n=y.n, w.n=w.n, cols=cols)
  sfStop()
  names(res) <- chrs
  res
}

consolidate.lows <- function(dat, cols=2:3, limit=0) {
  dat[rowSums(dat[,cols]) > limit,]
}

poiComparePerWin <- function(y, n.reads,grp=factor(rep(1:3, each=3))) {
 
  grp1 <- grp2 <- grp
  contrasts(grp1) <- cbind(c(0,1,0),
                           c(0,0,1))

  contrasts(grp2) <- cbind(c(0,1,1),
                           c(0,-1,1))
  dat1 <- data.frame(y=unlist(y[-1]), grp=grp1)
  dat2 <- data.frame(y=unlist(y[-1]), grp=grp2)

  ctl.idx <- which(grp==1)
  tst.ct1 <- glm(y~grp, data=dat1, family=quasipoisson,
                 offset=log(n.reads)-log(mean(n.reads[ctl.idx])))
  tst.ct2 <- glm(y~grp, data=dat2, family=quasipoisson,
                 offset=log(n.reads)-log(mean(n.reads[ctl.idx])))
  res1 <- as.vector(summary(tst.ct1)$coefficients)
  res2 <- as.vector(summary(tst.ct2)$coefficients[-1,])
  disp <- as.vector(summary(tst.ct1)$dispersion)
  overallP <- drop1(tst.ct1, test="F")$'Pr(F)'[2]
  res <- c(res1[1:3],res2[1:2], # estimate
           res1[4:6], res2[3:4], # std
           res1[7:9], res2[5:6], # t value
           res1[10:12], res2[7:8], # Pr
           disp, # dispersion
           overallP) # F test on grp
  res
}

par.poiComparePerWin <- function(dat.chr,
                                 dir,
                                 file.pre="",
                                 n.reads,
                                 grp=factor(rep(1:3, each=3)),
                                 n=10000,
                                 n.clusters=2,
                                 sourcefile=NULL)
{
  idx.split <- split(seq(dim(dat.chr)[1]), (seq(dim(dat.chr)[1])-1)%/%n)
  sfInit(parallel = T, cpus = n.clusters)
  sfExport('dat.chr', local=T)
  ##sfSource(sourcefile)
  names(idx.split) <- paste("s", names(idx.split),sep="")
  sfLapply(idx.split,
           function(idx, dat, n.reads, grp, dir, file.pre){
             filenm <- paste(dir,
                             paste(file.pre, paste(range(dat[idx,1]),collapse="-"),"RData", sep="."),
                             sep="/")
             if(file.access(filenm)<0) {
               res <- apply(dat[idx,], 1, poiComparePerWin, n.reads);
               res <- as.data.frame(t(res));
               row.names(res) <- dat[idx,1];
               save(res, file=filenm);
             }
           },
           dat=dat.chr, n.reads=n.reads, grp=grp, dir=dir, file.pre=file.pre)
  sfStop()
  ## res <- as.data.frame(t(as.data.frame(res)))
  ##   row.names(res) <- dat.chr[,1]
  ##   names(res) <- c(sapply(c("est", "std", "t", "P"),
  ##                        function(x, vec) paste(x, vec, sep="."),
  ##                        vec=c("int", "T1vsN","T2vsN","T1vsT2", "T1T2vsN")), "disp", "p.overall")
  ##   res
}

dir.poiComparePerWin <- function(dir,
                                 n.reads,
                                 grp=factor(rep(1:3, each=3)),
                                 n=100000,
                                 n.clusters=2) {
  files <- grep("JoinedSamples.chr", list.files("../data/Colon"), v=T)
  for (eachFile in files) {
    datnm <- load(paste(dir, eachFile, sep="/"))
    dat <- consolidate.lows(get(datnm), cols=2:10)
    if(dim(dat)[1]>0) {
      res <- par.poiComparePerWin(dat,
                                  dir=dir,
                                  file.pre=datnm,
                                  n.reads=n.reads,
                                  n=n,
                                  n.clusters=n.clusters)
      ## save(res, file=paste(dir, paste("poi",datnm, "RData", sep="."), sep="/"))
      ## cat(paste(datnm,"is done","\n"))
    }
  }
}

  


###. mixPoisson mcmc
##____________________________________________________________________

## FUNCTION fits 3-mixture poisson model with Gibbs algorithm
## 			lambda2-gamma(a,b), lambda3i-exp(r), r-gamma(s,q)
## ARGUMENT
##   y: count data
##   initials: initial values of model parameters if provided
##   controls: MCMC operating parameters
## OUTPUT
##   p.mcmc, lambda.mcmc, r: MCMC chain of model parameters
##   w.(win, var): membership probability (estimate, variance)
##   signal.(win, var): lambda3i (estimate, variance)
##   paraCRM: common estimates of p, lambda and r

mixPois_eWin <- function(y, initials=NULL, para.priors, controls,
                         p.tolerance=1e4)
{

  n <- length(y)
  if (min(y)<0|| any(is.na(y))) stop("Invalid count data y.")

  msize <- controls$msize
  burnin <- controls$burnin
  thin <- controls$thin

  d0 <- para.priors$d
  lambda.prior <- para.priors$lambda.prior
  a0 <- lambda.prior$a
  b0 <- lambda.prior$b
  tmp <- para.priors$signal.prior
  s0 <- tmp$s
  q0 <- tmp$q
  if (min(c(d0, a0, b0, q0, s0))<0) stop("Invalid prior distribution.")
  
  if (!is.null(initials)){
    p.old <- initials$p	
    lambda.old <- initials$lambda
    signal.old <- initials$signal
    r.old <- initials$r 
  } else {
    ## yn <- y[which(y>0)]; ycut <- quantile(yn, probs=0.5);
    ## p.old <- c(sum(y==0)/n, sum(y>0&y<=4)/n, sum(y>4)/n);
    ## lambda.old <- mean(y[which(y>0&y<=4)], na.rm=TRUE);
    ## r.old <- 1/mean(y[which(y>4)], na.rm=TRUE);
    ## signal.old <- y*(1+r.old);
    ###YXU: yn <- y[which(y>0)]; ycut <- quantile(yn, probs=0.5);
    yn <- y[y>0]
    ycut <- quantile(yn, probs=0.5)
    ###YXU: p.old <- c(sum(y==0)/n, sum(y>0&y<=4)/n, sum(y>4)/n);
    p.old <- c(sum(y==0)/n, sum(((y>0)+(y<=4))==2)/n, sum(y>4)/n)
    ###YXU: lambda.old <- mean(y[which(y>0&y<=4)], na.rm=TRUE);
    lambda.old <- mean(y[((y>0)+(y<=4))==2], na.rm=TRUE)
    ###YXU:r.old <- 1/mean(y[which(y>4)], na.rm=TRUE);
    r.old <- 1/mean(y[((y>0)+(y<=4))==2], na.rm=TRUE)
    signal.old <- y*(1+r.old)
  }
  ##browser()
  if (sum(p.old)-1>p.tolerance|min(p.old)<0|max(p.old)>1) stop("Invalid p initials.")
  if (min(lambda.old, r.old, signal.old)<0) stop("Invalid lambda initials.")
    
  p.mcmc <- matrix(NA, ncol=3, nrow=msize)
  lambda.mcmc <- rep(NA, msize) ## rate of measurement error noise;
  r.mcmc <- rep(NA, msize)
  
  w.pst <- matrix(0, ncol=3, nrow=n)
  signal.pst <- rep(0, n) ## posterior estimate of lambda3i
  w.var <- matrix(0, ncol=3, nrow=n)
  signal.var <- rep(0, n)
  
  ## mcmc iterations
  for (m in 1:msize){
    denpois <- cbind(dpois(y, lambda=0),
                     dpois(y, lambda=lambda.old),
                     dpois(y, lambda=signal.old))
    tmp <- denpois*p.old	
    w <- tmp/rowSums(tmp) 
    if (m > burnin){
      w.pst <- w.pst + w
      w.var <- w.var + w*w
    }
    ##browser()
    ## step 1 (Z)
    ##z.new <- apply(w, 1, rmultin) ## QQ
    ##z.new <- t(z.new)
    z.new <- matrix(0, nrow=n, ncol=3)
    uu <- runif(n, 0, 1)
    w2 <- w[,1]+w[,2]
    ###YXU: m123 <- ifelse(uu<=w[,1], 1, 2)
    m123 <- (uu>w[,1])+1
    ###YXU: tmp <- which(m123==2)
    tmpidx <- which(m123==2)
    ###YXU: m123[tmp] <- ifelse(uu[tmp]>w2[tmp], 3, 2)
    m123[tmpidx] <- ((uu>w2)+2)[tmpidx]
    ###YXU:
    z.new[which(m123==1),1] <- 1
    z.new[which(m123==2),2] <- 1
    z.new[which(m123==3),3] <- 1
    ### Not working in R12.1.: z.new <- cbind(m123==1,m123==2,m123==3)+0
    ## step 2 (P)
    dn <- d0 + colSums(w, na.rm=TRUE) 
    p.new <- rdirichlet(1, dn)	# row-vector
    ## step 3 (lambda)
    szy <- sum(z.new[,2]*y)
    sz <- sum(z.new[,2])
    lambda.new <- rgamma(1, shape=a0+szy, rate=b0+sz)
    ## step 4 (signal)
    szy <- z.new[,3]*y
    sz <- z.new[,3]
    ##signal.new <- rgam.main(shape=1+szy, rate=r.old+sz, lam=lambda.new)
    signal.new <- rgam.main2(shape=1+szy, rate=r.old+sz)
    ##signal.new <- rep(NA, n)	
    ##for (i in 1:n){ ## QQ	
    ##	tmp <- rgamma(1, shape=1+szy[i], rate=r.old+sz[i])
    ##	while(tmp<=lambda.new){ ## orderd 0<lambda<signali
    ##		tmp <- rgamma(1, shape=1+szy[i], rate=r.old+sz[i])
    ##		}
    ##	signal.new[i] <- tmp
    ##	}
    if (m>burnin){
      signal.pst <- signal.pst + signal.new
      signal.var <- signal.var + signal.new*signal.new
    }
    ## step 5 (r)
    hb1 <- (n+q0)+1
    hb2 <- sum(signal.new)+s0
    r.new <- rgamma(1, shape=hb1, rate=hb2)
    ## done
    p.mcmc[m,] <- p.new
    lambda.mcmc[m] <- lambda.new
    r.mcmc[m] <- r.new
    p.old <- as.vector(p.new)
    lambda.old <- lambda.new
    r.old <- r.new 
    signal.old <- signal.new

    if (m%%1000==0) print(paste("MCMC runs", m, "iterations."))
  }
  
  ## Posterior estimate of W (two methods)
  w.win <- w.pst/(msize-burnin)
  signal.win <- signal.pst/(msize-burnin)
  
  w.var <- w.var/(msize-burnin)-w.win*w.win
  signal.var <- signal.var/(msize-burnin)-signal.win*signal.win
  
  lambda.med <- median(lambda.mcmc)
  p.med <- apply(p.mcmc, 2, median)
  r.med <- median(r.mcmc)
  
  id <- seq(burnin, msize, by=thin)
  p.mcmc <- p.mcmc[id, ]
  lambda.mcmc <- lambda.mcmc[id]
  r.mcmc <- r.mcmc[id]
  
  return(list(w.win=w.win,
              signal.win=signal.win,
              w.var=w.var,
              signal.var=signal.var,
              paraCRM=list(p=p.med, lambda=lambda.med, r=r.med),
              chains=list(p.mcmc=p.mcmc, lambda.mcmc=lambda.mcmc, r.mcmc=r.mcmc)))
}


## Aux Functions 
rmultin <- function(p){	
  u <- runif(1, 0, 1)
  
  if (u<=p[1]){
    z <- c(1, 0, 0)
  } else {
    if (u<=(p[1]+p[2])){
      z <- c(0,1,0)	
    } else {
      z <- c(0,0,1)
    }
  }
  z
}




## DO NOT REQUIRE lambad3i>lambda2

rgam2 <- function(shape, rate){
  
  n <- length(shape)
  if (min(shape)<1) stop("Invalid shape in rgam().")
  
  d <- shape-1/3  
  c <- 1/(3*sqrt(d))
  
  x <- rnorm(n, 0, 1)
  u <- log(runif(n, 0, 1))
  v <- (1+c*x)*(1+c*x)*(1+c*x)
  
  ###YXU: accept <- ifelse(v>0, 1, 0)
  accept <- (v>0)+0
  tmp <- v>0
  cc <- 0.5*x[tmp]*x[tmp]+d[tmp]-d[tmp]*v[tmp] + d[tmp]*log(v[tmp])
  ###YXU: accept[tmp] <- ifelse(u[tmp]<=cc, 1, 0)
  accept[tmp] <- (u[tmp]<=cc)+0
  
  yc <- d*v
  y <- yc
  y[accept!=1] <- NA
  ##browser()
  ##y <- ifelse(accept==1, yc, NA)
  y <- y/rate
  ##y <- ifelse(y>lam, y, NA)
  y
}


rgam.main2 <- function(shape, rate){
  y <- rgam2(shape, rate)  
  ## Fill in NAs
  id <- is.na(y)
  this.shape <- shape[id]
  this.rate <- rate[id]
  y[id] <- rgam2(this.shape, this.rate)		
  
  ## Fill in NAs
  id <- is.na(y)
  this.shape <- matrix(shape[id], ncol=1)
  this.rate <- rate[id]

  y[id] <- rgamma(length(this.shape), shape=this.shape, rate=this.rate)
  ## for (i in 1:length(id)){
  ##   this.i <- id[i]
  ##   tmp <- rgamma(1, shape=this.shape[i], rate=this.rate[i])
  ##   y[this.i] <- tmp
  ## }	
  y
}

###. pattern recognition win
##____________________________________________________________________

###.. binomial
##--------------------------------------
b.f <- function(k, m, p) {
  dbinom(k, m, p)
}
Fb.f <- function(k,m,p) {
  res <- rep(NA, length(k))
  res[k<0] <- 0
  res[k>=0] <- pbinom(k[k>=0], m, p)
  res
}
  
Q2 <- function(k, m, p) {
  ## (4.4)
  Fb.f(k-1, m, p)^2
  - (k-1)*b.f(k, m, p)*Fb.f(k-2, m, p)
  + m*p*b.f(k, m, p)*Fb.f(k-3, m-1, p)
}

Q3 <- function(k, m, p) {
  stopifnot(k>2)
  A1 <- 2*b.f(k, m, p)*Fb.f(k-1, m, p) * ((k-1)*Fb.f(k-2,m,p)-m*p*Fb.f(k-3,m-1,p))
  A2 <- .5*b.f(k,m,p)^2*((k-1)*(k-2)*Fb.f(k-3,m,p)-2*(k-2)*m*p*Fb.f(k-4,m-1,p)+m*(m-1)*p^2*Fb.f(k-5,m-2,p))
  A3 <- sum(b.f(2*k-seq(k-1),m,p)*Fb.f(seq(k-1)-1,m,p)^2)
  A4 <- sum(b.f(2*k-seq(k-1),m,p)*b.f(seq(k-1),m,p)*((seq(k-1)-1)*Fb.f(seq(k-1)-2,m,p)-m*p*Fb.f(seq(k-1)-3,m-1,p)))
  Fb.f(k-1,m,p)^3-A1+A2+A3-A4
}

QL <- function(k,m,L,p) {
  ## N=m*L
  exp(log(Q2(k,m,p))+(L-2)*(log(Q3(k,m,p))-log(Q2(k,m,p))))
}

Pr.f <- function(k,m,L,p){
  1-QL(k,m,L,p)
}


f.f <- function(r,k,p) {
  p/r*dbinom(r-1,k-1,p)*(r*(1-p)*dbinom(r-1,k-1,p)+(r-k*p)*pbinom(r-2,k-1,p))
}

Pr.S.up <- function(r,k,n,p) {
  lambda <- (n-k+1)*f.f(r,k,p)
  exp(-1*lambda)
}


lambda <- function(r,k,n,p) {
  res <- (n-k+1)*(r/k-p)*dbinom(r,k,p)
}

Pr.S <- function(r,k,n,p) {
  exp(-lambda(r,k,n,p))
}



H.f <- function(theta, p) {
  theta*log(theta/p)+(1-theta)*log((1-theta)/(1-p))
}

h.f <- function(theta, p) {
  log(theta*(1-p)/(p*(1-theta)))
}

Pr.S.cp <- function(r,k,n,p, rho=1) {
  exp(-1*n*(r-k*p)*exp(-1*(k*H.f(r/k,p)+rho*h.f(r/k,p)))/sqrt(2*pi*r*k*(k-r)))
}


###.. Construct the cutoff table
##--------------------------------------

compute.cutoff <- function(r,n,p, cutoff=0.05) {
  k.seq <- seq(r,floor(r/p))
  k.idx <- which((1-Pr.S.up(r, k.seq, n, p) < cutoff))
  res <- k.seq[rev(k.idx)[1]]
  if(length(res)>0) res else NA
}
  
###.. Tree based search
##--------------------------------------
node.get.p.n <- function(anode, win.size) {
  res <- list(n=NA, p=NA)
  if(!is.leaf(anode)) {
    ewins <- as.numeric(labels(anode))/win.size
    n <- (max(ewins)-min(ewins))
    r <- attr(anode, 'members')
    p <- r/n
    res <- list(n=n, p=p)
  }
  res
}


node.compute.Pr <- function(anode, n, p, k.cutoff=100,win.size=100) {
  ##browser()
  if(!is.leaf(anode)) {
    ewins <- as.numeric(labels(anode))/win.size
    k <- (max(ewins)-min(ewins))
    r <- attr(anode, 'members')
    ##cat(r,'\t',k,'\t',n,'\t',p,'\n')
    if(k>=r/p | k>k.cutoff) 1 else 1-Pr.S.up(r,k,n,p)}
  else NA
}

node.set.Pr.aux <- function(anode, n, p, k.cutoff=100,win.size=100) {
  if(!is.leaf(anode)) {
    value <- node.compute.Pr(anode, n, p, k.cutoff, win.size) }
  else value <- NA
  attr(anode, 'Pr') <- value
  anode
}

node.set.Pr <- function(anode,n,p,k.cutoff=100, win.size=100) {
  dendrapply(anode, node.set.Pr.aux, n=n, p=p, k.cutoff=k.cutoff, win.size=win.size)
}


node.set.significance.aux <- function(anode, Pr.cutoff=0.05) {
  if(!is.leaf(anode)) {
    this.Pr <- attr(anode, 'Pr')
    ##if(attr(anode, 'Pr') < Pr.cutoff) attr(anode, 'Sig') <- TRUE else attr(anode, 'Sig') <- FALSE
    if(!is.na(this.Pr) && this.Pr < Pr.cutoff) attr(anode, 'Sig') <- TRUE
    else attr(anode, 'Sig') <- FALSE
  }
  else attr(anode, 'Sig') <- NA
  anode
}

node.set.significance <- function(anode, Pr.cutoff=0.001) {
  dendrapply(anode, node.set.significance.aux,Pr.cutoff=Pr.cutoff)
}

node.pick.range <- function(anode, id.start=1,tol=1e4, debug=F){
  rests <- NULL
  i <- id.start
  getmembers <- function(X) {
    if(!is.leaf(X)) {
      output <- F
      pr0 <- attr(X, 'Pr')
      pr1 <- attr(X[[1]], 'Pr')
      pr2 <- attr(X[[2]], 'Pr')
      lf1 <- is.leaf(X[[1]])
      lf2 <- is.leaf(X[[2]])
      if(debug) cat(pr0,pr1,pr2,lf1,lf2)
      if(attr(X, 'Sig')) {
        if(lf1 & !lf2 & pr2 > pr0) output <- T
        else if(!lf1 & lf2 & pr1 > pr0) output <- T
        else if(lf1 & lf2) output <- T
        else if(!lf1 & !lf2 & (pr0-min(pr1, pr2))<tol) output <- T
      }

      if(debug) cat('-->',output,'\n')

      if(output) {
        this.wins <- sort(labels(X))
        this.ID <- rep(i, length(this.wins))
        this.pr <- rep(attr(X, "Pr"), length(this.wins))
        this.res <- cbind(this.ID, this.wins, this.pr)
        rests <<- rbind(rests,this.res)
        i <<- i+1
        
      }
      else {
        for (jj in seq(length=length(X))) {
          X[[jj]] <- getmembers(X[[jj]]) 
        }
      }
    }
    X
  }
  getmembers(anode)
  if(!is.null(rests)) {
    rests <- as.data.frame(rests)
    names(rests) <- c('ID','Win','Pr')
  }
  rests
}

events.score <- function(res, 
                         events.res, 
                         rownames, 
                         sample.label, 
                         which.score=2, 
                         npart=10, 
                         para.mode=2) {
  ## so called event probability

  if(length(which.score)>1) stop("Select only one number from 1-5!")
  events.res <- as.data.frame(events.res)
  sel.nonull <- !unlist(lapply(res, is.null)); res <- res[sel.nonull]
  
  p3 <- c(unlist(lapply(res, get.p)))
  w3.res <- as.matrix(t(laply(res, get.w3)))
  dimnames(w3.res) <- NULL
  cat("num of input sites:", nrow(w3.res),"\n")
  e.res <- vector('list', ncol(events.res))
 
  patts <- pattern.design(events.res, sample.label)
  for (i in seq(ncol(events.res))) {
    cat("pattern:",paste(patts[,i],sep="",collapse=""))
    .this <- event.prob(w3.res, 
                        sample.label,
                        event=patts[,i],
                        p=p3, 
                        which.score=which.score, 
                        npart=npart,
                        para.mode=para.mode)
    e.res[[i]] <- as.vector(.this)
    cat("\t","done\n")
  }
  e.res <- as.data.frame(e.res)
  names(e.res) <- paste("P", apply(events.res,2,paste, collapse=""), sep="")
  row.names(e.res) <- rownames
  attr(e.res, 'p3') <- p3
  attr(e.res, 'event') <- events.res
  attr(e.res, 'label') <- sample.label
  e.res
}

events.cutoff <- function(res, e.res, cutoff, method = c("beta", "BF", "best"),find.range= c(-20,20)) {
  ## Beta prior based distribution
  ## events =
  ##   0 0
  ##   1 0
  ##   1 1
  ## events are defined by columns
  ##browser()

  method <- match.arg(method)
  events <- as.data.frame(attr(e.res, 'event'))
  sample.label <- attr(e.res, 'label')
  if(method=="beta") {
    res <- res[!unlist(lapply(res, is.null))]
    para.names <- paste("P",apply(events, 2, paste, collapse=""), sep="")
    events <- as.data.frame(pattern.design(events, sample.label))
    para <-  apply(events,
                   2,
                   function(x) compute.minB.para(res, x)
                   )
    names(para) <- para.names
    
    e.cut <- lapply(para,
                    function(x) find.cutoff(cutoff,
                                            x,
                                            find.range)
                    )
    
    ret <- unlist(e.cut)
    attr(ret, 'para') <- para
    names(ret) <- NULL
  }
  else if(method=="BF") {
    ret <- rep(cutoff, ncol(events))
    attr(ret, 'para') <- NA
  }
  else if(method=="best") {
    .cutoff <- findBest.escoreCutoff(e.res, seq(-2,2,0.01), cutoff)
    ret <- rep(.cutoff, ncol(events))
    attr(ret, 'para') <- NA
  }
  attr(ret, 'method') <- method
  ret
}

choose.events.cutoff <- function(res, e.res, cutoff, method=c("BF","beta","best")) {
  method <- match.arg(method)
  ## if(method == "BF") out <- events.cutoff(res, e.res, cutoff, 'BF')
  ## else if(method == "beta") out <- events.cutoff(res, e.res, cutoff, 'beta')
  out <- events.cutoff(res, e.res, cutoff, method)
  out
}


pattern.design <- function(pattern, sample.label) {
  sl <- factor(sample.label)
  patt <- as.data.frame(pattern)
  stopifnot(nlevels(sl) == nrow(patt))
  levels(sl) <- seq(nlevels(sl))
  ret <- apply(patt, 2, function(x) x[sl])
  if(is.null(dim(ret))) ret <- matrix(ret, nrow=1)
  dimnames(ret) <- NULL
  ret
}

which.zeroPattern <- function(pattNames, cap="P") {
  which(as.numeric(gsub(cap, "", pattNames))==0)
}

win.select <- function(e.res, e.cut, exclude.zeroPattern=T, zeroPattern.cutoff=0){
  ##browser()
  method = attr(e.cut, 'method')
  stopifnot(length(e.res)==length(e.cut))
  ## select the windows by thresholding on event probabilities or BF
  p3 <- attr(e.res, 'p3')
  para <- attr(e.cut, "para")
  events.res <- attr(e.res, 'event')
  n.win <- nrow(e.res)
  n.event <- ncol(e.res)
  wins <- row.names(e.res)

  wins.sel <- vector('list', n.event)
  for (i in seq(n.event)) {
    sel.i <- e.res[,i]> e.cut[i]
    scores.i <- e.res[sel.i,i]
    wins.i <- as.numeric(wins[sel.i])
    p.i <- vector('numeric', length(scores.i))
    if(method=="beta" & length(scores.i)>0) p.i <- sapply(scores.i,compute.minB.p, para[[i]])
    else {
      ecdf.i <- ecdf(e.res[,i])
      p.i <- 1-ecdf.i(scores.i)
    }
    wins.sel[[i]] <- data.frame(scores=scores.i, wins=wins.i, probs=p.i,row.names=wins.i)    
  }
  names(wins.sel) <- names(e.res)
  attr(wins.sel, 'event') <- events.res
  attr(wins.sel, 'cutoff') <- e.cut
  ## removing zero patterns from other if asked
  if(exclude.zeroPattern) {
    idx.0 <- which.zeroPattern(names(e.res))
    idx.oth <- setdiff(seq(n.event), idx.0)
    wins.zero <- wins[e.res[,idx.0]>zeroPattern.cutoff]
    for (j in idx.oth) {
      wins.sel[[j]] <- subset(wins.sel[[j]],subset=!(wins.sel[[j]]$wins%in%wins.zero))
    }
  }
  wins.sel
}

win.overlap.pair <- function(wins.sel.pair) {
  wins.sel.1 <- wins.sel.pair[[1]]
  wins.sel.2 <- wins.sel.pair[[2]]
  win.int <- intersect(wins.sel.1$wins, wins.sel.2$wins)
  del.in.2 <- wins.sel.2[win.int,"probs"] >= wins.sel.1[win.int,"probs"]
}

identify.win.overlap <- function(wins.sel) {
  
  ow <- NULL
  seq.n <- seq(length(wins.sel))
  for (i in seq.n ) {
    ow <- c(ow, wins.sel[[i]]$wins)
  }

  o.wins <-  names(which(table(ow)>1))

  ow.p <- NULL
  for (i in seq.n) {
    ow.p <- cbind(ow.p, wins.sel[[i]][o.wins,"probs"])
  }

  zero.col <- which(colSums(attr(wins.sel, 'event')) == 0)
  if(length(zero.col) ==1) {
    ow.p.this <- ow.p[,-zero.col,drop=F]
    if(ncol(ow.p.this)>1) {
      ids <- apply(ow.p.this,1, which.min)
    } else {
      ids <- NULL
    }
  }
  else if (length(zero.col) ==0) {
    if(ncol(ow.p)>1) {
      ids <- apply(ow.p,1, which.min)
    } else {
      ids <- NULL
    }
  }
  else stop("Wrong number of columns for zero patterns!")
  ##ids.full <- apply(ow.p,1, which.min)
  ret <- vector("list", length(wins.sel))
  for (i in seq.n) {
    o.wins.i <- o.wins[ids != i]
    wins.sel.i <- wins.sel[[i]]
    ret[[i]] <- wins.sel.i[as.character(setdiff(wins.sel.i$wins, o.wins.i)),]
  }
  attributes(ret) <- attributes(wins.sel)
  ret
}
  
get.wiSum <- function(res) {
  ret <- lapply(res, function(x) if(!is.null(x)) colSums(x$w.win) else NULL)
  ret <-  t(as.data.frame(ret[!unlist(lapply(ret, is.null))]))
  row.names(ret) <- NULL
  ret
}

get.w3Chr <- function(res) {
  ret <- lapply(res, function(x) if(!is.null(x)) x$paraCRM$p[3] else NULL)
  c(unlist(ret))
}

compute.beta.para <- function(res, d = c(1/3, 1/3, 1/3)) {
  wisum <- get.wiSum(res)
  d <- matrix(rep(d, each=nrow(wisum)), ncol= length(d))
  wd <- wisum+d
  cbind(wd[,3], wd[,1]+wd[,2])
}

compute.beta.para.2 <- function(res) {
  a <- get.w3Chr(res)
  cbind(a=a, b=1-a)
}

compute.beta.para.3 <- function(res, rev=F) {
  if(rev) {
    ret <- lapply(res, function(x) if(!is.null(x)) est.beta.para(1-x$w.win[,3]))
  }
  else {
    ret <- lapply(res, function(x) if(!is.null(x)) est.beta.para(x$w.win[,3]))
  }
  ret <- as.data.frame(matrix(unlist(ret), byrow=T, ncol=2))
  names(ret) <- c("a", "b")
  ret
}

compute.minB.para <- function(res,
                              pattern=c(0,0,0,1,1,1,1,1,1),
                              v=3)
{
  if(v==1) ab.1 <- compute.beta.para(res)
  else if(v==2) ab.1 <- compute.beta.para.2(res)
  else  ab.1 <- compute.beta.para.3(res)
  e.1 <- get.w3Chr(res)
  if(length(pattern) != nrow(ab.1) | length(pattern) != length(e.1)) stop("dimension mismatch!")
  ab.0 <- ab.1[,rev(seq(ncol(ab.1)))]
  e.0 <- 1 - e.1
  
  ab <- ab.1
  e <- e.1

  ab[pattern == 1,] <- ab.1[pattern == 1,]
  ab[pattern == 0,] <- ab.0[pattern == 0,]

  e[pattern == 1] <- e.1[pattern == 1]
  e[pattern == 0] <- e.0[pattern == 0]

  list(ab=ab, e=e)
}

compute.minB.p <- function(x, paras) {
  a <- paras$ab[,1]
  b <- paras$ab[,2]
  e <- paras$e
  q <- (2^x) * e/(1 - e + (2^x)*e)
  pb <- pbeta(q, a, b)
  ##print(pb)
  prod(1-pb)
}

est.beta.para <- function(x) {
  ## Moment estimate of parameters
  u <- mean(x)
  v <- var(x)
  a <- (1-u)*u^2/v-u
  b <- (1-u)*(u*(1-u)/v-1)
  list(a=a, b=b)
}

find.cutoff.objective <- function(x, paras, cutoff) {
  (compute.minB.p(x=x, paras=paras)-cutoff)^2
}

find.cutoff <- function(cutoff, paras, range=c(-10,10)) {
  ret <- optimize(find.cutoff.objective,
                  range,
                  paras=paras,
                  cutoff=cutoff)
  round(ret$minimum, 3)
}
         
  


win.den <- function(wins.sel.i, win.size=100) {

  #sel.wins <- as.numeric(wins[e>e.cut])
  #sel.wins.df <- data.frame(sel.wins, row.names=sel.wins)
  sel.wins.df <- wins.sel.i[,2,drop=F]
  den.res <- as.dendrogram(hclust(dist(sel.wins.df/win.size), method='complete'))
  den.res.reorder <- reorder(den.res, sel.wins.df$wins, mean)
  den.res.reorder
}


win.search <- function(atree, cut=4e5, win.size=100) {
  atree.l <- cut(atree,cut)$lower
  res.wins <- NULL
  id.count <- 1
  allWins <- labels(atree)
  for (i in seq(atree.l)) {
    np <- node.get.p.n(atree.l[[i]],win.size)
    ##if (i >1 & nrow(res.wins) >0) id.count <- max(as.numeric(res.wins$ID))+1
    acut <- node.pick.range(node.set.significance(node.set.Pr(atree.l[[i]],
                                                              np$n, np$p,
                                                              win.size=win.size)),
                            id.start=(if(i==1) id.count else max(as.numeric(res.wins$ID))+1))
                            ##id.start = id.count)
    res.wins <- rbind(res.wins, acut)
  }
  res.wins$Win <- as.character(res.wins$Win)
  res.wins$ID <- as.integer(res.wins$ID)
  res.wins$Pr <- as.numeric(as.character(res.wins$Pr))
  singleWins <- setdiff(allWins, res.wins$Win)
  singleWins.res <- data.frame(ID=seq(singleWins)+length(unique(res.wins$ID)),
                               Win=as.character(singleWins), Pr=rep(0, length(singleWins)),
                               stringsAsFactors=F)
  rbind(res.wins, singleWins.res)
}

match.mx <- function(x.mx, table.mx) {
  if(is.null(dim(table.mx))) dim(table.mx) <- c(1,length(table.mx))
  reslt <- rep(0, nrow(x.mx))
  reslt <- unlist(apply(x.mx, 1, function(a.x.mx, table.mx) {
                            rst <- which(apply(table.mx, 1, 
                                        function (a.table.mx, a.x.mx) all(a.x.mx == a.table.mx),
                                        a.x.mx = a.x.mx                           
                                       )
                                );
                            rst <- ifelse(length(rst)>0, rst,NA);
                            rst
                            },
                table.mx = table.mx
               ),
               
               )
  reslt
}


search.sites <- function(wins.sel, events.test, chr, win.size=100, den.cut=5e5){
  sites <- NULL
  for.seq <- match.mx(t(events.test), t(attr(wins.sel, 'event')))
  site.names <- names(wins.sel)
  for (i in for.seq) {
    den.i <- win.den(wins.sel[[i]],
                     win.size=win.size)
    sites.i <- win.search(den.i, cut=den.cut, win.size=win.size)
    chr.i <- rep(chr, nrow(sites.i))
    patt.i <- rep(site.names[i], nrow(sites.i))
    
    sites <- rbind(sites, cbind(chr.i, patt.i, sites.i))
  }
  names(sites) <- c("chr", "pattern", names(sites.i))
  attr(sites, 'event') <- events.test
  sites
}

search.sites.v2 <- function(wins.sel,
                            events.test=NULL,
                            chr,
                            win.size=100,
                            den.cutoff = c(1e5, 5e4),
                            cut.1st=c(3000,2000))
{
### omit empty wins.sel
  noemp <- sapply(wins.sel, function(x) nrow(x)>0)
  attr.event <- attr(wins.sel, 'event')[,noemp,drop=F]
  attr.cutoff <- attr(wins.sel, 'cutoff')[noemp]
  wins.sel <- wins.sel[noemp]
  attr(wins.sel, 'event') <- attr.event
  attr(wins.sel, 'cutoff') <- attr.cutoff
 
  if(is.null(events.test)) for.seq <- seq(wins.sel)
  else {
    wins.sel <- pattern.events(wins.sel, events.test)
    for.seq <- seq(wins.sel)
  }
  
  den.cutoff <- den.cutoff[attr(wins.sel,'nonull'),,drop=F]
  events.test <- events.test[,attr(wins.sel,'nonull'),drop=F]
  ##for.seq <- match.mx(t(events.test), t(attr(wins.sel, 'event')))

  site.names <- names(wins.sel)
  if(is.null(nrow(den.cutoff))) den.cutoff <- matrix(rep(den.cutoff,
                                                         each=length(for.seq)),
                                                     ncol=2)
  else if(nrow(den.cutoff) != length(for.seq)) stop("cutoff dim does not match event dim!")

  sfExport("wins.sel","events.test","chr","win.size","den.cutoff","cut.1st")
  wins.sel.cuts <- sfSapply(for.seq,
                            function(x) cut.wins(wins.sel[[x]], den.cutoff[x,1],cut.1st[1],cut.1st[2]),
                            simplify=F)
  wins.sel.cuts.size <- sapply(wins.sel.cuts,
                               seq,
                               simplify=F)
  idx <- NULL
  for(i in seq(along=wins.sel.cuts)) {
    idx <- rbind(idx, cbind(i, seq(length(wins.sel.cuts[[i]]))))
  }
  idx.list <- as.list(as.data.frame(t(idx)))
  process.search <- function(l) {
    ##browser()
    den.cut.2 <- den.cutoff[l[1],2]
    wins.sel.cut <- wins.sel.cuts[[l[1]]][[l[2]]]
    if(nrow(wins.sel.cut)>1) {
      den.i <- win.den(wins.sel.cut,win.size=win.size)
      sites.i <- win.search(den.i, cut=den.cut.2,win.size=win.size)
    }
    else if(nrow(wins.sel.cut) == 1) {
      sites.i <- data.frame(ID=1, Win=wins.sel.cut$wins, Pr=0)
    }
    else stop('nrow(wins.sel.cut) == 0')
    sites.i$ID <- paste(l[2], sites.i$ID, sep=".")
    chr.i <- rep(chr, nrow(sites.i))
    patt.i <- rep(site.names[l[1]], nrow(sites.i))
    cbind(chr.i, patt.i, sites.i)
  }
  ##sfSource(src.file)
  
  out.list <- sfClusterApplyLB(idx.list, process.search)
  out <- NULL
  for (i in seq(along=out.list)) {
    out <- rbind(out, out.list[[i]])
  }
  names(out) <- c("chr", "pattern", names(out)[-c(1,2)])
  attr(out, 'patterns') <- events.test
  attach.winScoreP(out, wins.sel)
}

attach.winScoreP <- function(sites, wins.sel) {
  out <- NULL
  for ( ipatt in names(wins.sel)) {
    patt.sub <- subset(sites, subset=pattern==ipatt)
    patt.score.p <- wins.sel[[ipatt]][patt.sub$Win,]
    out <- rbind(out, patt.score.p)
  }
  outall <- cbind(sites, out[,c('scores','probs')])
  names(outall) <- c(names(sites), c('Win.score', 'Win.P'))
  outall
}
    


cut.wins <- function(wins.df, h, cutoff.cat=3000, cutoff.last=2000)
{
  wins <- wins.df$wins
  tst <- diff(wins)
  ##browser()

  wins.cuts <- unique(c(0, wins[which(tst > h)], max(wins)+1))
  wins.cat <- cut(wins, wins.cuts, right=F)

  tb <- as.data.frame(table(as.numeric(wins.cat)))

  id <- rep(0, nrow(tb))
  cnt <- 0
  restart <- F
  for (i in seq(nrow(tb)))  {
    if(i==1) {
      id[i] = tb$Var1[i]
      cnt = tb$Freq[i] }
    else {
      if(!restart) {
        cnt <- cnt + tb$Freq[i]
        id[i] <- id[i-1]
      }
      else {
        cnt <- tb$Freq[i]
        id[i] <- max(id)+1
      }
      
      if(cnt > cutoff.cat) restart <- T
      else restart <- F
    }
  }
  is.one.group <- length(unique(id)) == 1
  if(cnt < cutoff.last & !is.one.group) id[id==id[which.max(id)]] <- max(id)-1
  
  tabl <- cbind(tb, ID=id)
  wins.cat.2 <- as.numeric(wins.cat)
  for (i.ID in unique(id)) {
    wins.cat.2[as.numeric(wins.cat) %in% tb[id==i.ID,'Var1']] <- i.ID
  }

  ##browser()
  id.unique <- unique(id)
  res <- vector('list', length=length(id.unique))
  names(res) <- id.unique
  for (i in seq(along=id.unique)) {
    sel <- wins.cat.2==id.unique[i]
    if(sum(sel) >0) res[[i]] <- wins.df[sel,]
  }
  res    
}

cut.wins2 <- function(wins.df, h, winsize) {
  ### partially done what's in cut.wins with better performance
  ### smaller trees are not merged into the bigger branch
  ### rewrite win.search
  ### NOT READY FOR USE
  .ir <- IRanges(start=wins.df$wins, width=winsize)
  .ir.gaps <- gaps(.ir)
  .gaps.cut <- .ir.gaps[width(.ir.gaps)>h]
  split.fac <- cut(wins.df$wins,
                   breaks=sort(c(min(wins.df$wins),
                     start(.gaps.cut),
                     end(.gaps.cut),
                     max(wins.df$wins))),
                   include.lowest=T)
  split(wins.df, split.fac, drop=T)
}

get.reads <- function(reads, start, end) {
  idx.s <- as.numeric(reads$poi) >= as.numeric(start)
  idx.e <- as.numeric(reads$poi) <= as.numeric(end)
  reads[idx.s & idx.e,]
}


get.reads.poi <- function(reads.poi, start, end) {
  idx.s <- as.numeric(reads.poi) > as.numeric(start)
  idx.e <- as.numeric(reads.poi) < as.numeric(end)
  reads.poi[idx.s & idx.e]
}

shear.sites <- function(sites, gap, reads.poi, wins.sel, win.size=100) {
  stopifnot(gap>win.size)
  gap <- gap/win.size
  for (i in unique(sites$pattern)) {
    sites.i <- subset(sites, subset=pattern==i)
    ids.i <- unique(sites.i$ID)
    ## KEEP IT HERE
    ##zwins <- c(wins.sel[[as.character(i)]][,'wins'], wins.sel[['P000']][,'wins'])
    ##zwins <- c(as.integer(sites.i[,'Win']), wins.sel[['P000']][,'wins'])
    ##zwins <- as.integer(sites.i[,'Win'])
    ## KEEP IT HERE
    for ( j in ids.i) {
      sites.i.j <- subset(sites.i, subset=ID==j)
      ids.tb <- table(sites.i.j$ID)
      ids <- names(ids.tb[ids.tb>1])
      sites.ijo <- sites.i.j[order(as.numeric(sites.i.j$Win)),]
      sites.ijo.wins <- as.integer(sites.ijo$Win)
      low.ijod <- sites.ijo.wins[-1]
      up.ijod <- sites.ijo.wins[-nrow(sites.ijo)]
      sites.ijod <- data.frame(up=up.ijod,
                               low=low.ijod,
                               d=((low.ijod-up.ijod)/win.size)> gap)
      sites.ijod.sel <- subset(sites.ijod, subset=d)
      ## KEEP IT HERE before confirmation
      ## gap.TF <- apply(sites.ijod.sel,
      ##                 1,
      ##                 function(x, zw) !all(get.reads.poi(reads.poi, x[1], x[2])%in% zw), zw=zwins)
      ## KEEP IT HERE
      gap.TF <- rep(TRUE,nrow(sites.ijod.sel))
      
      sites.ijod.sel.gap <- cbind(sites.ijod.sel, gap=gap.TF)
      if(nrow(sites.ijod.sel.gap) > 0) {
        ids.ijo <- sites.ijo$ID
        ids.ijo.out <- ids.ijo
        wins.ijo <- as.integer(sites.ijo$Win)
        aid <- 1
        
        for ( k in seq(nrow(sites.ijod.sel.gap))) {
          
          gap.i <- sites.ijod.sel.gap[k,'gap']
          gap.w <- sites.ijod.sel.gap[k,'up']
          if(gap.i) {
            ids.ijo.out[wins.ijo > gap.w] <- paste(ids.ijo[wins.ijo > gap.w], aid, sep=".")
            aid <- aid + 1
          }
        }

        sites.ijo$ID <- ids.ijo.out
        sites[sites$pattern==i & sites$ID ==j,] <- sites.ijo
      }
    }
  }
  sites
}

join.sites <- function(sites, gap=500, win.size=100) {
  ##TODO: a c version
  stopifnot(gap>win.size)
  gap <- gap/win.size
  patts <- unique(sites$pattern)
  res <- sites
  for (i in patts) {
    sub.i <- which(sites$pattern == i)
    sites.i <- sites[sub.i,]
    sites.by <- by(as.numeric(sites.i$Win),sites.i$ID, range)

    sites.by.i <- as.data.frame(matrix(unlist(sites.by), ncol=2, byrow=T))
    names(sites.by.i) <- c('start','end')
    row.names(sites.by.i) <- dimnames(sites.by)[[1]]

    sites.by.i.o <- sites.by.i[order(sites.by.i[,1]),]
    sites.by.i.o.IDs <- row.names(sites.by.i.o)
    sites.by.i.d <- data.frame(lw=sites.by.i.o.IDs[-1],
                               up=sites.by.i.o.IDs[-nrow(sites.by.i)],
                               d=((sites.by.i.o[-1,1]-sites.by.i.o[-nrow(sites.by.i),2])/win.size)<= gap)
    sites.by.i.d.sel <- sites.by.i.d[which(sites.by.i.d$d),]

    if (nrow(sites.by.i.d.sel) > 0) {
      tst.id <- as.integer(row.names(sites.by.i.d.sel))
      tst.lw <- as.character(sites.by.i.d.sel$lw)
      tst.up <- as.character(sites.by.i.d.sel$up)

      if(nrow(sites.by.i.d.sel) > 1) {
        for ( i in 2:length(tst.id)) {
          if((tst.id[i]-tst.id[i-1])==1) tst.up[i] <- tst.up[i-1]
        }
      }
    
      for ( j in seq(length(tst.lw))) {
        sites.i$ID[sites.i$ID==tst.lw[j]] <- tst.up[j]
      }
    }
    
    res[sub.i,] <- sites.i
  }
  res
}

## Define patterns that are of interest for searching

pattern.events <- function(wins.sel, events.comp) {
  ## events.comp
  # NA NA
  # NA 1
  # 1  1
  ## use NA for elements that are free of change
  ##browser()
  events.comp <- as.matrix(events.comp)
  events <- attr(wins.sel, "event")
  comp.list <- NULL
  out.list <- NULL
  for (i in seq(ncol(events.comp))) {
    rows2sel <- !is.na(events.comp[,i])
    events.i <- as.matrix(events[rows2sel,,drop=F])
    events.sel <- which(!is.na(match.mx(t(events.i),
                                        t(events.comp[rows2sel,i,drop=F])
                                        )))
    comp.list <- c(comp.list, list(events.sel))
  }
  for (j in seq(comp.list)) {
    this.comp <- NULL
    for(l in comp.list[[j]]) {
      this.comp <- rbind(this.comp, wins.sel[[l]])
    }
    out.list <- c(out.list, list(this.comp))
  }
  out.list.names <- paste("P",apply(events.comp, 2, paste, collapse=""), sep="")
  out.list.names <- gsub("NA", "-", out.list.names)
  names(out.list) <- out.list.names
  ### take care the null list
  nonullidx <- !sapply(out.list, is.null)
  out.list <- out.list[nonullidx]
  ###
  attr(out.list, 'event') <- events.comp[,nonullidx,drop=F]
  attr(out.list, 'nonull') <- nonullidx
  out.list
}

match.pattern.win <- function(wins.sel, sites) {
  ##browser()
  pattern.win <- rep("", nrow(sites))
  sites.Win <- sites$Win
  wins.sel.names <- names(wins.sel)
  for(i in seq(along=wins.sel)) {
    idx <- sites.Win %in% wins.sel[[i]]$wins
    ## idx.noNA <- na.omit(idx)
    ## if(length(idx.noNA) > 0) {
    ##   print(length(idx.noNA))
    pattern.win[idx] <- paste(pattern.win[idx], wins.sel.names[i], sep="")
    ## }
  }
  pattern.win
}

size.site <- function(sites) {
  ## browser()
  output <- NULL
  sites.unique.pattern <- unique(sites$pattern)
  for (i in seq(along=sites.unique.pattern)) {
    this.sites <- subset(sites, subset=pattern==sites.unique.pattern[i])
    sites.Win <- as.numeric(this.sites$Win)
    sites.ID <- as.character(this.sites$ID)
    ##sites.pattern <- sites$pattern
    this.output <- tapply(sites.Win,
                          sites.ID,
                          function(x) max(x)-min(x),
                          simplify=T)
    idx <- match(unique(sites.ID), names(this.output))
    output <- rbind(output, data.frame(ID=unique(sites.ID),
                                       size=this.output[idx],
                                       pattern=sites.unique.pattern[i],
                                       stringsAsFactors=F))
  }
  
  output
}

event.allPossiblePatterns <- function(smplabel) {
  alabel <- smplabel
  grp.n <- length(unique(alabel))
  as.data.frame(t(expand.grid(as.data.frame(matrix(rep(c(0,1),grp.n), ncol=grp.n)))))
}

event.patternsToScan <- function(allpatt, exclude=NULL,noZeroPatt=T) {
  
  out <- allpatt
  if(!is.null(exclude)) out <- out[,-1*exclude,drop=F]
  if(noZeroPatt & length(which(colSums(out) == 0))>0) out <- out[,-which(colSums(out) == 0),drop=F]
  out
}

# ###. Initial filtering
# ##____________________________________________________________________
# est.pois.lambda <- function(tab) {
#   f1 <- with(tab, freq[count==1])
#   f2 <- with(tab, freq[count==2])
#   2*f2/f1
# }

# est.fdr <- function(n,tb) {
#   lambda <- est.pois.lambda(tb)
#   n.total <- with(tb, freq[count==1])/dpois(1, lambda)
#   n.total*dpois(n,lambda)/with(tb, freq[count == n])
# }

# compute.initfilter.cutoff <- function(tbs, cutoff=0.5) {
#   ## which count cutoff should be used so that at least one sample start to contain at least 'cutoff'
#   ##  positives (max of counts across all samples should be at east this number)
#   ifsatisfy <- rowSums(sapply(tbs, function(x) sapply(1:10, est.fdr, tb=x)) < cutoff) >0
#   which(ifsatisfy)[1]
# }

# get.initfilter.cutoff <- function(counts.tab, cutoff=0.5) {
#   ##get.mappedCountTab()
#   compute.initfilter.cutoff(counts.tab, cutoff=cutoff)
# }
###. Normalization
##____________________________________________________________________
est.NB <- function(tb) {
  lambda <- est.pois.lambda(tb)
  with(tb, freq[count==1])/dpois(1, lambda)
}

## norm.ratio <- function(tb, lambda.ref) {
##   Nb <- est.NB(tb)
##   lambda <- est.pois.lambda(tb)
##   cnt <- tb$count
##   1+Nb*(dpois(cnt, lambda.ref)- dpois(cnt, lambda))/tb$freq
## }

## norm.freq <- function(tb, lambda.ref) {
##   r <- norm.ratio(tb, lambda.ref)
##   tb$freq <- tb$freq*r
##   tb
## }

bg.N <- function(tbs) {
  bg.lambda <- mean(sapply(tbs, est.pois.lambda), trim=.1)
  bg.N <- median(sapply(tbs, est.NB))
  max.N <- max(sapply(tbs, function(x) max(x$count)))
  count <- seq(0, max.N)
  freq <- bg.N*dpois(count, bg.lambda)
  data.frame(count, freq)
}
  

denoiseSignal <- function(tb) {
  ##browser()
  N <- sapply(tb, est.NB)
  lams <- sapply(tb, est.pois.lambda)
  newtb <- vector("list", length(tb))
  for (i in seq(along=tb)) {
    newtb[[i]] <- tb[[i]]
    newtb[[i]]$freq <- as.integer(tb[[i]]$freq-N[i]*dpois(tb[[i]]$count, lams[i]))
  }
  newtb
}

norm.tbs <- function(tbs, bg.correct=F) {
  ##browser()
  tb.dn <- denoiseSignal(tbs)
  bg <- bg.N(tbs)
  ## tb.norm <- lapply(tb.dn, function(x) {x$freq <- x$freq + bg[x$count+1,'freq'];
  ##                                       x=rbind(c(0, bg[bg$count==0, 'freq']),x);
  ##                                     x})
  if(bg.correct) {
    tb.norm <- tb.dn
  } else {
    tb.norm <- lapply(tb.dn, function(x) {x$freq <- x$freq + bg[x$count+1,'freq'];x})
  }
   
  tb.norm
}

extend.tbs <- function(tbs) {
  lapply(tbs, function(x) rbind(c(0,est.NB(x)*dpois(0, est.pois.lambda(x))), x))  
}
             
normalize.signal <- function(x, tb, tb.norm, out.int=F) {
  ##browser()
  ls0 <- approxfun(tb$count, tb$freq)
  lsnorm <- approxfun(tb.norm$count,tb.norm$freq)

  out <- x*lsnorm(x)/ls0(x)
  if(out.int) out <- as.integer(round(out, 0))
  out
}

normalize.pmeans <- function(count, countab, bg.correct=F, out.int=F) {
  ##browser()
  cnt.ex <- extend.tbs(countab)
  cnt.norm <- norm.tbs(cnt.ex, bg.correct=bg.correct)
  norm.sig <- as.data.frame(sapply(seq(ncol(count)), function(i) normalize.signal(count[,i],
                                                                                  cnt.ex[[i]],
                                                                                  cnt.norm[[i]],
                                                                                  out.int=out.int)))
  row.names(norm.sig) <- row.names(count)
  names(norm.sig) <- names(count)
  norm.sig
}

subset.pmeans <- function(pmeans, subset) {
  list(mean=pmeans$mean[,subset],
       var=pmeans$var[,subset])
}
                        

###. Differential Testing
##____________________________________________________________________

get.sig <- function(x) {
  x$signal.win
}

get.sigvar <- function(x) {
  x$signal.var
}

get.w <- function(x, which=3) {
  x$w.win[,which]
}

get.p <- function(x, which=3) {
  x$paraCRM$p[which]
}

get.w.var <- function(x) {
  x$w.var[,3]
}

compute.p <- function(x, event=c(0,0,0,1,1,1,1,1,1))
  ## assume x is in w3.res format as above
{
  res <- x
  for (i in seq(event)) {
    if(event[i] !=1) res[,i] <-  1-x[,i]
  }
  res
}

compute.v <- function(p,p.var,k=1e-36) {
  var.zeros <- min(p.var[p.var !=0])
  p.var[p.var ==0] <- var.zeros
  ##p.zeros <- min(p[p!=0])
  ##p.one <- max(p[p!=1])
  ((p+k)^2*(1-p+k)^2)/((1-p+k)^2+(p+k)^2)*1/p.var
}

compute.T <- function(p) {
  logit.f(p)
}

compute.Tbar <- function(T,v) {
  res <- T
  for (i in seq(ncol(T))) {
    res[,i] <- T[,i]*v[,i]
  }
  rowSums(res)/rowSums(v)
}

get.postMean <- function(x, log=F) {
  if(log) {
    res <- log(x$signal.win+0)*x$w.win[,3]+log(x$paraCRM$lambda+0)*x$w.win[,2]
  }
  else res <- x$signal.win*x$w.win[,3]+x$paraCRM$lambda*x$w.win[,2]
  res
}

get.postLogLambda2 <- function(x, log=F) {
  if(log) {
    res <- log(x$paraCRM$lambda+0)*x$w.win[,2]
  }
  else res <- x$paraCRM$lambda*x$w.win[,2]
}

get.postLogLambda22 <- function(x, log=F) {
  if(log) {
    res <- log(x$paraCRM$lambda+0)*x$w.win[,2]^2
  }
  else res <- x$paraCRM$lambda*x$w.win[,2]^2
}

get.lambdaVar <- function(x) {
  if(is.null(x)) res <- NULL
  else res <- var(x$chains$lambda.mcmc)
  res
}

get.postLogLambda3 <- function(x, log=F) {
  if(log) {
    res <- log(x$signal.win+0)*x$w.win[,3]
  }
  else res <- x$signal.win*x$w.win[,3]
}

get.postLambda3.var <- function(x) {
  x$signal.var*(x$w.win[,3])^2
}

get.postmean.var <- function(x) {
  x$signal.var/(x$signal.win)^2
}

get.lambda <- function(x) {
  x$paraCRM$lambda
}

get.postMean.cutoffs.aux <- function(pmeans, lambda, cut.p=0.05) {
  pm1 <- pmeans[pmeans<lambda]
  pm2 <- pmeans[pmeans>lambda]

  cut1 <- quantile(pm1, 1-cut.p)
  cut2 <- quantile(pm2, cut.p)
  c(cut2, cut1)
}

get.postMean.cutoffs.aux2 <- function(pmeans, lambda, cut.p=0.05) {
  pm1 <- pmeans[pmeans<lambda]
  pm2 <- pmeans[pmeans>lambda]

  cut1 <- quantile(pm1, 1-cut.p)
  cut2 <- quantile(pm2, cut.p)
  (cut1+cut2)/2
}

get.pmeans.cut <- function(pmeans.df,lambda.vec, cut.p=0.05, method=2) {
  stopifnot(nrow(pmeans.df)!=length(lambda.vec))
  res <- NULL
  for ( i in seq(lambda.vec)) {
    if(method!=2) {
      res <- rbind(res, get.postMean.cutoffs.aux(pmeans.df[,i], lambda.vec[i], cut.p=cut.p))
    }
    else {
      this.res <- get.postMean.cutoffs.aux2(pmeans.df[,i], lambda.vec[i], cut.p=cut.p)
      res <- rbind(res, c(this.res, this.res))
    }
  }
  res <- as.data.frame(res)
  names(res) <- c('lo', 'hi')
  res
}
    

assign.groups <- function(pmeans, cuts, patt=c(0,0,0,1,1,1,1,1,1)){
  res <- pmeans
  for (i in seq(ncol(pmeans))) {
    res[,i] <- if(patt[i] ==0) pmeans[,i]<= cuts$lo[i] else pmeans[,i] > cuts$hi[i]
  }
  apply(res,1, all)
}

get.siteSig <- function(sigres, msites, which.pattern, which.site) {
  Patt.idx <- grep(which.pattern,msites$PatternID)
  Win.idx <- which(msites$SiteID==which.site)
  wins <- msites[intersect(Patt.idx, Win.idx),'Win']
  sigres[as.character(wins),]
}

logit.f <- function(x, k=1e-36) {
  log((x+k)/(1-x+k))
}

get.w3 <- function(x) {
  x$w.win[,3]
}

get.w3.var <- function(x) {
  x$w.var[,3]
}

event.prob.mx <- function(x, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), k=1e-36)
  ## assume x is in w3.res format as above
  ## p is the p3
{
    
  res <- x
  for (i in seq(along=event)) {
    if(event[i] !=1) {
      x[,i] <-  1-x[,i]
      p[i] <- 1-p[i]
    }
  }

  for (i in seq(along=event)) {
    res[,i] <- log2((x[,i]+k)/(1-x[,i]+k)*((1-p[i])/p[i])) ## test on p + k
  }
  res
}

## event.prob <- function(x, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), k=1e-36, quantile.p=0.25, quantile.type=7)
##   ## assume x is in w3.res format as above
##   ## p is the p3
## {
##   res <- event.prob.mx(x, event, p, k)
##   ## res <- x
##   ## for (i in seq(along=event)) {
##   ##   if(event[i] !=1) {
##   ##     x[,i] <-  1-x[,i]
##   ##     p[i] <- 1-p[i]
##   ##   }
##   ## }

##   ## for (i in seq(along=event)) {
##   ##   res[,i] <- log2((x[,i]+k)/(1-x[,i]+k)*((1-p[i])/p[i])) ## test on p + k
##   ## }  
  
##   #output <- apply(res, 1, min)
##   output <- apply(res, 1, function(x) quantile(x, probs=quantile.p, type=quantile.type, names=F,na.rm=T))
##   output
## }
event.prob <- function(x, sample.label,event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), k=1e-36,
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

findBest.escoreCutoff <- function(escores, cuts=seq(-2,2,0.01), cutoff=0.1) {
  n.tot <- nrow(escores)[1]
  miss.n <- sapply(cuts, function(k, dat) sum(rowSums(dat > k) > 1), dat=escores)
  miss.pct <- miss.n/n.tot
  cuts[which(miss.pct<=cutoff)[1]]
}

escore.comp <- function(x, group, which.p) {
  min(tapply(x,
             group,
             function(x, which.p) quantile(x, p=which.p),
             which.p=which.p))
}

event.prob.cut <- function(n, event=c(0,0,0,1,1,1,1,1,1), p=rep(0.33,9), p.cut=0.001,k=1e-36) {
  res <- NULL
  for (i in seq(event)) {
    if(event[i] !=1) {
      p[i] <- 1-p[i]
    }
    res <- cbind(res, rbeta(n, p[i], 1-p[i]))
  }
  
  

  for (i in seq(event)) {
    res[,i] <- log2((res[,i]+k)/(1-res[,i]+k)*((1-p[i])/p[i])) ## test on p + k
  }
  quantile(apply(res, 1, min), 1-p.cut)
}
  
event.v.prob <- function(x,v,  event=c(0,0,0,1,1,1,1,1,1)) {
  p <- x
  for (i in seq(event)) {
    if(event[i] !=1) p[,i] <-  1-x[,i]
  }
  logp <- log(p)
  logw <- v/p^2
  logres <- logp*logw
  res <- apply(logres,1,sum)
  res
}

format.data.lme <- function(dat, group=c(1,1,1,2,2,2,3,3,3), smpid=NULL, log.t=TRUE) {
  y <- as.numeric(unlist(dat))+0
  if(log.t) y <- log2(y)
  win <- rep(row.names(dat),ncol(dat))
  smp <- rep(names(dat), each=nrow(dat))
  grp <- as.factor(rep(group, each=nrow(dat)))
  res <- data.frame(y=y, win=win, smp=smp, grp=grp)
  if(!is.null(smpid)) {
    pid <- as.factor(rep(smpid, each=nrow(dat)))
    res <- cbind(res, smpid=smpid)
  } 
  res
}

## test.effect <- function(wins, psig, vars, contr= cbind(c(0,1,1),c(0,1,-1)),
##                         multicore=T, n.core=4)
test.effect <- function(wins, psig, vars, group.label,
                        n.core,
                        contr= cbind(c(0,1,1),c(0,1,-1)),
                        src.file="")
{
  ##browser()
  unique.wins <- unique(wins$ID)
  n.wins <- length(unique.wins)
  test.effect.worker <- function(i) {
    ##browser()
    ## 31775 is one with NULL output
    ### Debug code
    ## if(i%%10 ==0) cat("** ")
    ## cat(i)
    ## if(i%%10==0) cat("\n")

    wins.by.id <- as.character(unlist(subset(wins,
                                             subset=ID==unique.wins[i],
                                             select="Win")))
    pmean.sub <- psig[wins.by.id,]
    pmean.sub.lme <- format.data.lme(pmean.sub,
                                      group=as.numeric(factor(group.label)),
                                      log.t=TRUE)
    pmean.var.sub <- vars[wins.by.id,]
    pmean.sub.lme.wt <- data.frame(pmean.sub.lme,
                                   w = as.vector(unlist(pmean.var.sub)))
      
    contrasts(pmean.sub.lme.wt$grp) <- make.RContrast(contr)  
    varfun <- varIdent(form = ~ 1 | grp)
    lme.pmean.sub.wt <- try(lme(y ~ grp,
                                data=pmean.sub.lme.wt,
                                random= ~ 1 |win,
                                control=list(
                                  msMaxIter=500,
                                  maxIter=500,
                                  opt="optim",
                                  tolerance=1e-3),
                                weights=varfun),
                            TRUE)
    if(!is(lme.pmean.sub.wt, 'try-error')) {
      sigmaSq <- diag((lme.pmean.sub.wt$sigma*coef(lme.pmean.sub.wt$modelStruct$varStruct,
                                                   uncons=F,
                                                   allCoef=TRUE))^2)
      ##R00R00 <- lme.pmean.sub.wt$varFix%*%ginv(sigmaSq)
      R00R00 <- lme.pmean.sub.wt$varFix%*%diag(1/diag(sigmaSq))
      
      sigmaSq.d <- diag(as.vector(by(pmean.sub.lme.wt$w, pmean.sub.lme.wt$grp, mean)))
      
      sigmaSq.updated <- sigmaSq + sigmaSq.d
      
      varFix.updated <- R00R00%*%sigmaSq.updated
      
      stdErr.updated <- sqrt(diag(varFix.updated))
      ## browser()
      this.tTable <-  summary(lme.pmean.sub.wt)$tTable
      this.tTable[,"Std.Error"] <- stdErr.updated
      this.tTable[,"t-value"] <- this.tTable[,"Value"]/stdErr.updated
      this.tTable[,"p-value"] <- pt(-abs(this.tTable[,"t-value"]), this.tTable[,"DF"])*2
      res.i <- this.tTable
      retn <- as.vector(res.i)
      nm.1 <- dimnames(res.i)[[1]]
      nm.2 <- dimnames(res.i)[[2]]
      retn.names <- as.vector(t(sapply(nm.1, function(x) paste(x, nm.2, sep="."))))
      names(retn) <- retn.names      
    } else {
      retn <- NULL
    }
    retn
  }
  if(n.core>1) n.wins.cut <- as.numeric(cut(seq(n.wins), n.core)) else n.wins.cut <- rep(1, n.wins)
  n.wins.seq <- seq(n.wins)

  ##sfExportAll()
  ##sfSource(src.file)
  sfExport("wins", "psig", "vars", "contr", "n.core", "n.wins.cut","n.wins.seq",
           "unique.wins","n.wins")

  ##sfSource()
  test.effect.worker.group <- function(i) {
    ## cat(i,"\n")
    res.i <- sapply(n.wins.seq[n.wins.cut==i],test.effect.worker, simplify=F)

    ## debug:
    ## res.i <- sapply(n.wins.seq[n.wins.cut==i],
    ##                 function(x) {trythis <- try(test.effect.worker(x))
    ##                              if(is(trythis, "try-error")) cat("*",x, "\n")
    ##                              trythis
    ##                            },
    ##                 simplify=F)
    ##t(res.i)
    res.i
  }
  ##res <- unlist(sfLapply(seq(n.core), test.effect.worker.group), recursive=F)
  res.list <- sfLapply(seq(n.core), test.effect.worker.group)

  ## take the 1st non-null for names and ncols
  res.names <- NULL
  res.ncols <- NULL
  for(i in seq(along=res.list)){
    for(j in seq(along=res.list[[i]])) {
      a.res <- res.list[[i]][[j]]
      if(!is.null(a.res)) {
        res.names <- names(a.res)
        res.ncols <- length(a.res)
        break
      }    
    }
    if(!is.null(res.ncols)) break
  }

  ##browser()
  res <- as.data.frame(matrix(unlist(res.list), ncol=res.ncols, byrow=T))
  names(res) <- res.names
  
  ##res.notNull <- !as.vector(unlist(lapply(res.list, is.null)))
  res.notNull <- !as.vector(unlist(lapply(res.list, function(x) sapply(x, is.null))))
  ##res.list <- res.list[res.notNull]
  row.names(res) <- unique.wins[res.notNull]
  res <- data.frame(ID=as.character(unique.wins[res.notNull]), res, stringsAsFactors=F)
  res
}

test.effect.rep <- function(wins, psig, vars, group.label, rep.label,
                            n.core,
                            contr= cbind(c(0,1,1),c(0,1,-1)),
                            src.file="")
{
  ##browser()
  unique.wins <- unique(wins$ID)
  n.wins <- length(unique.wins)
  test.effect.worker <- function(i) {
    ##browser()
    ## 31775 is one with NULL output
    ### Debug code
    ## if(i%%10 ==0) cat("** ")
    ## cat(i)
    ## if(i%%10==0) cat("\n")
    wins.by.id <- as.character(unlist(subset(wins,
                                             subset=ID==unique.wins[i],
                                             select="Win")))
    pmean.sub <- psig[wins.by.id,]
    pmean.sub.lme <- format.data.lme(pmean.sub,
                                     group=as.numeric(factor(group.label)),
                                     smpid=rep.label,
                                     log.t=TRUE)
    pmean.var.sub <- vars[wins.by.id,]
    pmean.sub.lme.wt <- data.frame(pmean.sub.lme,
                                   w = as.vector(unlist(pmean.var.sub)))
    contrasts(pmean.sub.lme.wt$grp) <- make.RContrast(contr)  
    varfun <- varIdent(form = ~ 1 | grp)
    lme.pmean.sub.wt <- try(lme(y ~ grp,
                                data=pmean.sub.lme.wt,
                                random= ~ 1 |smpid/win,
                                control=list(
                                  msMaxIter=500,
                                  maxIter=500,
                                  opt="optim",
                                  tolerance=1e-3),
                                weights=varfun),
                            TRUE)
    if(!is(lme.pmean.sub.wt, 'try-error')) {
      sigmaSq <- diag((lme.pmean.sub.wt$sigma*coef(lme.pmean.sub.wt$modelStruct$varStruct,
                                                   uncons=F,
                                                   allCoef=TRUE))^2)
      ##R00R00 <- lme.pmean.sub.wt$varFix%*%ginv(sigmaSq)
      R00R00 <- lme.pmean.sub.wt$varFix%*%diag(1/diag(sigmaSq))

      sigmaSq.d <- diag(as.vector(by(pmean.sub.lme.wt$w, pmean.sub.lme.wt$grp, mean)))
      
      sigmaSq.updated <- sigmaSq + sigmaSq.d
      
      varFix.updated <- R00R00%*%sigmaSq.updated
      
      stdErr.updated <- sqrt(diag(varFix.updated))
      ## browser()
      this.tTable <-  summary(lme.pmean.sub.wt)$tTable
      this.tTable[,"Std.Error"] <- stdErr.updated
      this.tTable[,"t-value"] <- this.tTable[,"Value"]/stdErr.updated
      this.tTable[,"p-value"] <- pt(-abs(this.tTable[,"t-value"]), this.tTable[,"DF"])*2
      res.i <- this.tTable
      retn <- as.vector(res.i)
      nm.1 <- dimnames(res.i)[[1]]
      nm.2 <- dimnames(res.i)[[2]]
      retn.names <- as.vector(t(sapply(nm.1, function(x) paste(x, nm.2, sep="."))))
      names(retn) <- retn.names      
    } else {
      retn <- NULL
    }
    retn
  }
  if(n.core>1) n.wins.cut <- as.numeric(cut(seq(n.wins), n.core)) else n.wins.cut <- rep(1, n.wins)
  n.wins.seq <- seq(n.wins)
  ##sfExportAll(except=c( "globalNoExport" ))
  sfExport(list=ls())
  ##sfSource(src.file)
  ## sfExport("wins", "psig", "vars", "contr", "n.core", "n.wins.cut","n.wins.seq",
  ##          "unique.wins","n.wins","test.effect.worker","group.label","rep.label")
  ##sfSource()
  test.effect.worker.group <- function(i) {
    cat(i,"\n")
    res.i <- sapply(n.wins.seq[n.wins.cut==i],test.effect.worker, simplify=F)

    ## debug:
    ## res.i <- sapply(n.wins.seq[n.wins.cut==i],
    ##                 function(x) {trythis <- try(test.effect.worker(x))
    ##                              if(is(trythis, "try-error")) cat("*",x, "\n")
    ##                              trythis
    ##                            },
    ##                 simplify=F)
    ##t(res.i)
    res.i
  }
  ##res <- unlist(sfLapply(seq(n.core), test.effect.worker.group), recursive=F)
  
  res.list <- sfLapply(seq(n.core), test.effect.worker.group)
  ##res.list <- sapply(seq(n.core), test.effect.worker.group)
  
  ## take the 1st non-null for names and ncols
  res.names <- NULL
  res.ncols <- NULL
  for(i in seq(along=res.list)){
    for(j in seq(along=res.list[[i]])) {
      a.res <- res.list[[i]][[j]]
      if(!is.null(a.res)) {
        res.names <- names(a.res)
        res.ncols <- length(a.res)
        break
      }    
    }
    if(!is.null(res.ncols)) break
  }

  ##browser()
  res <- as.data.frame(matrix(unlist(res.list), ncol=res.ncols, byrow=T))
  names(res) <- res.names
  
  ##res.notNull <- !as.vector(unlist(lapply(res.list, is.null)))
  res.notNull <- !as.vector(unlist(lapply(res.list, function(x) sapply(x, is.null))))
  ##res.list <- res.list[res.notNull]
  row.names(res) <- unique.wins[res.notNull]
  res <- data.frame(ID=as.character(unique.wins[res.notNull]), res, stringsAsFactors=F)
  res
}

process.diffTest <- function(sites.js.ex, events.wins, events.ctr, testmode.rep,
                             pmeans.norm, pmeans.var,
                             ncore,
                             sample.label,
                             subject.label) {
  patts <- unique(sites.js.ex$pattern)
  events.wins.matched <- events.wins[,match.patt.pattmx(patts, events.wins)]
  events.ctr.matched <- events.ctr[match.patt.pattmx(patts, events.wins)]
  patt.n <- length(patts)
  wins.test <- vector("list", length=patt.n)
  cat("testing ...\n")
  if(testmode.rep) {
    for (i in seq(patt.n)) {
      wins.test[[i]] <- list()
      for(j in seq(length(events.ctr.matched[[i]]))) {
        cat(as.character(patts[i]), "\t", j, "\n")
        wins.test[[i]][[j]] <- test.effect.rep(subset(sites.js.ex, subset=pattern==patts[i]),
                                               pmeans.norm,
                                               pmeans.var,
                                               contr= events.ctr.matched[[i]][[j]],
                                               n.core = ncore,
                                               group.label=sample.label,
                                               rep.label=subject.label)
        attr(wins.test[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
        attr(wins.test[[i]][[j]],'pattern') <- events.wins.matched[,i]
      }
    }
  } else {
    for (i in seq(patt.n)) {
      wins.test[[i]] <- list()
      for(j in seq(length(events.ctr.matched[[i]]))) {
        cat(as.character(patts[i]), "\t", j, "\n")
        wins.test[[i]][[j]] <- test.effect(subset(sites.js.ex, subset=pattern==patts[i]),
                                           pmeans.norm,
                                           pmeans.var,
                                           contr= events.ctr.matched[[i]][[j]],
                                           n.core = ncore, group.label=sample.label)
        attr(wins.test[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
        attr(wins.test[[i]][[j]],'pattern') <- events.wins.matched[,i]
      }
    }
  }

  cat("report:\n")
  wins.test.report <- vector("list", length=patt.n)
  for (i in seq(patt.n)) {
    wins.test.report[[i]] <- list()
    for(j in seq(length(events.ctr.matched[[i]]))) {
      cat(as.character(patts[i]), "\t", j, "\n")
      wins.test.report[[i]][[j]] <- combine.report(subset(sites.js.ex, subset=pattern==patts[i]),
                                                   wins.test[[i]][[j]])
      attr(wins.test.report[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
      attr(wins.test.report[[i]][[j]],'pattern') <- events.wins.matched[,i]
    }
  }
  list(wins.test=wins.test, wins.test.report=wins.test.report)
}

combine.report <- function(wins, tests, tests.select) {
  ##browser()
  if(missing(tests.select)) tests.select=grep("ID|Value|p.value", names(tests))
  part1 <- as.data.frame(t(as.data.frame(lapply(by(as.integer(wins$Win),wins$ID, range), unlist))))
  rownames.part1 <- gsub('X', '', row.names(part1))
  row.names(part1) <- rownames.part1

  tests.sel <- tests[,tests.select]
  common.rownames <- intersect(rownames.part1, row.names(tests.sel))
  part2 <- tests.sel[common.rownames,]
  
  res <- data.frame(part1[common.rownames,], part2)
  ##row.names(res) <- common.rownames
  names(res) <- c(c("start","end"), names(part2))
  res
}


normalize.sig <- function(sig, bg) {
  centr <- apply(bg, 2, function(x) mean(x))
  diff <- centr-mean(centr)
  scl <- apply(bg, 2, function(x) sd(x))
  rescl <- scl/mean(scl)
  ##browser()
  for(i in seq(ncol(sig))) {
    sig[,i] <- (sig[,i]+diff[i])/rescl[i]
  }
  sig
}

normalize.sig.2 <- function(sig, bg, lambda22) {

  a <- apply(bg, 2, function(x) var(x))/apply(lambda22, 2, function(x) mean(x))
  rescl.a <- a/mean(a)

  for(i in seq(ncol(bg))) {
    bg[,i] <- bg[,i]/rescl.a[i]
  }
  
  centr <- apply(bg, 2, function(x) median(x))
  diff <- centr-median(centr)
  
  for(i in seq(ncol(sig))) {
    sig[,i] <- sig[,i]/rescl.a[i]-diff[i]
  }
  sig
}


select.sites <- function(reports,
                         cut.fdr=0.05, cut.effect = 1,
                         cut.base=0.2, cut.eq=0.2,
                         col.select=NULL, win.size=100) {

  if(is.null(col.select)) col.select <- seq(ncol(reports))
  sel.fdr <- p.adjust(reports$pvalue, 'fdr') < cut.fdr
  sel.effect <- abs(reports$effect) > cut.effect
  sel.base <- reports$pbase > cut.base
  sel.eq <- reports$peq > cut.eq

  sel <- apply(cbind(sel.fdr, sel.effect, sel.base, sel.eq), 1, all)
  res <- reports[sel,col.select]
  res$end <- res$end+win.size
  size <- res$end-res$start
  res <- cbind(res[,c(1,2)],size,res[,-c(1,2)])
  reorder <- order(res$pvalue)
  res[reorder,]
}

selReads.report <- function(id, wins, reads) {
  if('poi' %in% names(reads)) {
    reads.poi <- reads$poi
    reads <- reads[,-1]
    row.names(reads) <- reads.poi
  }

  id <- as.character(id)
  idx <- unlist(sapply(id, function(x,y) which(x==y), y=wins$ID))
  win.sel <- wins[idx,c('ID', 'Win')]
  read.sel <- reads[as.character(win.sel$Win),]
  res <- cbind(win.sel, read.sel)
  res
}

est.grpMeans <- function(all.test.report) {
  out.list <- all.test.report
  for (i in seq(all.test.report)) {
    for (j in seq(all.test.report[[i]])) {
      this.df <- all.test.report[[i]][[j]]
      a <- this.df[,'beffect']
      b <- this.df[,'effect']
      c <- this.df[,'eqffect']
      this.ex <- data.frame(baseline=a - b - c,
                            grp1=a + 2 * c,
                            grp2=a + b - c)
      out.list[[i]][[j]] <- cbind(this.df, this.ex)
     }
  }
  out.list
}

default.pattContrast <- function(patt) {
  ##browser()
  if(!all(patt==1)) {
    n.patt <- length(patt)
    patt.flip <- rep(1,n.patt)-patt
    if(patt.flip[1]>0) base.contr <- patt
    else base.contr <- patt.flip
    ## if(n.patt>2)
    ##   out.contr <- cbind(base.contr,
    ##                      matrix(0, nrow=n.patt, ncol=n.patt-2))
    ## else out.contr <- as.matrix(base.contr)
    out.contr <- as.matrix(base.contr)
    dimnames(out.contr) <- NULL
   
  }
  else out.contr <- NULL
  out.contr
}

make.RContrast <- function(pattContr) {
  ##browser()
  base.contr <- t(pattContr)[1,]
  zero.idx <- base.contr==0
  one.idx <- base.contr==1
  base.contr[zero.idx]= -1/sum(zero.idx)
  base.contr[one.idx]=1/sum(one.idx)
  out.contr <- make.contrasts(base.contr)
  dimnames(out.contr) <- NULL
  out.contr
}

make.reportContrast <- function(report.flat,
                                pvalue.var="p.value.1",
                                effect.var="value.1",
                                type=c(1,2))
{
  if(type ==1) {
    contrDF <- expand.grid(dimnames(with(report.flat, table(contrast, pattern))))
    contrList <- vector("list", nrow(contrDF))
    for(i in seq(length(contrList))) {
      contrList[[i]]$contrast <- as.character(contrDF$contrast[i])
      contrList[[i]]$pattern <- as.character(contrDF$pattern[i])
      contrList[[i]]$pvalue.var <- pvalue.var
      contrList[[i]]$effect.var <- effect.var
    }
  }
  else if(type==2) {
    acontr <- dimnames(with(report.flat, table(contrast, pattern)))
    acontr$pvalue.var <- pvalue.var
    acontr$effect.var <- effect.var
    contrList <- list(acontr)
  }
  contrList
}
    


default.contrastList <- function(patts) {
  ##browser()
  inpatt <- as.matrix(patts)
  out <- vector("list",ncol(inpatt))
  for(i in seq(dim(inpatt)[2])) {
    out[[i]] <- default.pattContrast(inpatt[,i])
  }
  unique(out)
}

construct.contrasts <- function(patts) {
  ##browser()
  patt.ncol <- ncol(as.matrix(patts))
  patt.sum <- colSums(patts)
  out <- vector("list", patt.ncol)
  for(i in seq(patt.ncol)) {
    ipatt <- patts[,i]
    if(patt.sum[i] != nrow(patts)) out[[i]] <- default.contrastList(ipatt)
    else {
      out[[i]] <- default.contrastList(patts[,patt.sum!=nrow(patts)])
    }
  }
  out
}

construct.contrast.objs <- function(events.wins, events.ctr) {
  unique.ctrs <- unique(unlist(lapply(events.ctr,
                                      function(x) lapply(x,
                                                         function(x) paste(apply(as.matrix(x),
                                                                                 2,
                                                                                 paste,
                                                                                 collapse=""),
                                                                           collapse="_")))))
  patts.str <- apply(events.wins, 2, paste, collapse="")
  sapply(unique.ctrs, function(x) list(contrast=x,
                                       pattern=patts.str,
                                       pvalue.var="p.value.1",
                                       effect.var="value.1"),
         simplify=F)
}


    
    
  

###. process report
##____________________________________________________________________

###.. Functions
##--------------------------------------
process.report <- function(chr, winsize=100) {
  regionPatt.file <- paste(region.path,
                           paste(chr,"consolidated_wins.RData", sep="."),
                           sep="/")
  mixPoi.file <- paste(mixPoi.path,
                       paste("mixPoi",chr,"RData", sep="."),
                       sep="/")
  diff.file <- paste(diff.path,
                     paste("differential",chr,"RData", sep="."),
                     sep="/")
  load(regionPatt.file)
  load(mixPoi.file)
  load(diff.file)

  patts <- unique(sites.js$pattern)
  patt.n <- length(patts)

  cat("*** ",chr," ***","\n")
  sel.test.report <- NULL
  sel.read.report <- NULL
  for (i in seq(patt.n)[-7]) {
    for(j in seq(length(events.ctr[[i]]))) {
      sel.test.report.0 <- select.sites(all.test.report[[i]][[j]],1,1,0,0,win.size=winsize)
      if(nrow(sel.test.report.0) > 0) {
        sel.read.report.0 <- selReads.report(row.names(sel.test.report.0),
                                             subset(sites.js, subset=pattern==patts[i]),
                                             reads)
  
        sel.test.report.0 <- data.frame(chr=chr,
                                        pattern=as.character(patts)[i],
                                        contrast=paste("c",
                                          paste(events.ctr[[i]][[j]][,1],
                                                collapse=""),
                                          sep=""),
                                        ID=row.names(sel.test.report.0),
                                        sel.test.report.0)
        sel.read.report.0 <- data.frame(chr=chr,
                                        pattern=as.character(patts)[i],
                                        contrast=paste("c",
                                          paste(events.ctr[[i]][[j]][,1],
                                                collapse=""),
                                          sep=""),
                                        sel.read.report.0)
        
        sel.test.report <- rbind(sel.test.report, sel.test.report.0)
        
        sel.read.report <- rbind(sel.read.report, sel.read.report.0)
        cat(as.character(patts[i]),"\t",j,"\t", nrow(sel.test.report.0),"\n")
      }
    }
  }
  list(sites=sel.test.report, reads=sel.read.report)
}

###... Compute overlaps
if.overlap <- function(r1, r2) {
  ret <- T
  this.order <- order(c(as.numeric(r1), as.numeric(r2)))
  this.order.rev <- order(c(as.numeric(r2), as.numeric(r1)))
  if(sum(this.order-seq(4) != 0) ==0) ret <- F
  else if(sum(this.order.rev-seq(4) != 0) == 0) ret <- F
  ret
}

any.overlaps <- function(R1, R2, e=500) {
  R1 <- t(apply(R1, 1, `+`, c(-e, e)))
  R2 <- t(apply(R2, 1, `+`, c(-e, e)))
  ret <- NULL
  for ( i in seq(nrow(R1))) {
    ret <- c(ret, any(apply(R2, 1, function(y,x) if.overlap(x,y), x=R1[i,])))
  }
  ret
}

any.overlap.sites <- function(R1, R2, e=500, self=F) {

  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  R1.ex <- R1
  R1.ex[,1] <- R1[,1]-e
  R1.ex[,2] <- R1[,2]+e
  overlap.count <- rep(0, nrow(R1))
  overlap <- vector("logical", nrow(R1))
  overlap.where <- rep(NULL, nrow(R1))

  for (i in seq(nrow(R1))) {
    r1 <- R1.ex[i,]
    if(self) {
      overlap.v <- ((R2[-i,1]-r1[2]) * (R2[-i,2]-r1[1])) <0
      overlap.count[i] <- sum(overlap.v)
      overlap[i] <- overlap.count[i]>0
      overlap.where[i] <- paste(which(overlap.v), collapse=",")
    }
    else {
      overlap.v <- ((R2[,1] - r1[2]) * (R2[,2] - r1[1])) < 0
      overlap.count[i] <- sum(overlap.v)
      overlap[i] <- overlap.count[i]>0
      overlap.where[i] <- paste(which(overlap.v), collapse=",")
    }
    
  }
  data.frame(overlap=overlap, count=overlap.count, where=overlap.where)
}

any.overlap.sites.allchr <- function(R1,R2,
                                     e=500, self=F,
                                     chr.col="chr",
                                     start.col="start",
                                     end.col="end") {
  chrs <- unique(R1[,chr.col])
  overlapped.R1.sites <- NULL
  for (chr.i in chrs) {
    R1.sub <- subset(R1,
                     subset=R1[,chr.col]==chr.i,
                     select=c(start.col, end.col))
    R2.sub <- subset(R2, subset=R2[,chr.col]==chr.i,
                     select=c(start.col, end.col))
    olp <- any.overlap.sites(R1.sub,
                             R2.sub,
                             e=e,
                             self=self)
    row.names(olp) <- row.names(R1.sub)
    overlapped.R1.sites <- rbind(overlapped.R1.sites, olp)
    
  }
  overlapped.R1.sites
}


overlap.proportion <- function(R1,R2) {
  ## R1 and R2 are expected to be matrix with two columns and same number of rows.
  
  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  
  delta1 <- R1[,2]-R2[,1]
  delta2 <- R1[,1]-R2[,2]
  olp <- delta1*delta2<0
  deltap <- apply(cbind(delta1, delta2, R1[,2]-R1[,1], R2[,2]-R2[,1]), 1, 
                  function(x) min(abs(x)))
  p1 <- deltap/(R1[,2]-R1[,1])*olp
  p2 <- deltap/(R2[,2]-R2[,1])*olp
  data.frame(p1=p1,p2=p2)
}

overlap.size <- function(R1,R2) {
  ## R1 and R2 are expected to be matrix with two columns and same number of rows.
  
  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  
  delta1 <- R1[,2]-R2[,1]
  delta2 <- R1[,1]-R2[,2]
  olp <- delta1*delta2<0
  deltap <- apply(cbind(delta1, delta2, R1[,2]-R1[,1], R2[,2]-R2[,1]), 1, 
                  function(x) min(abs(x)))
  deltap
}




##overlap.sites <- 

gen.poi <- function(poi.data, sizes, win.size=100)
{ ## generate the positions
  ##browser()
  pois.start <- sample(poi.data,length(sizes))
  pois.end <- pois.start+sizes
  R1 <- cbind(pois.start, pois.end)
  Rex <- R1
  Rex[,1] <- R1[,1]-2000
  Rex[,2] <- R1[,2]+2000
  ##overlaps <- any.overlaps(R1, R1, e=2000)
  overlaps <- vector("logical", nrow(R1))
  for(i in seq(nrow(Rex))) {
    r1 <- Rex[i,]
    overlaps[i] <- any((Rex[-i,1]-r1[2]) * (Rex[-i,2]-r1[1]) <0)
  }
  R1 <- R1[!overlaps,]
  apply(R1, 1, function(x, by) seq.int(from=x[1], to=x[2], by=by), by=win.size)
}

sel.sigSites <- function(a.test.report, method="fdr", cutoff=0.05) {
  sigsites <- subset(a.test.report,
                     subset=p.adjust(a.test.report$pvalue, method) < cutoff)
  sig.pos <- subset(sigsites, effect > 0)
  sig.neg <- subset(sigsites, effect < 0)
  list(up=sig.pos, down=sig.neg)
}

###.. Multiple testing adjustment and site selection
##--------------------------------------
select.sites.GWadj <- function(sites.data,
                               reads.data,
                               method=c("qvalue","BY","bonferroni"),
                               cutoff=0.05,
                               pvalue.colname="pvalue",
                               sites.out.file="allsites",
                               reads.out.file="allreads",
                               out.path="Reports")
{
  method <- match.arg(method)
  out.sites.sel <- NULL
  out.reads.sel <- NULL
  if(method == "qvalue") {
    q.sites <- compute.qvalues(sites.data[,pvalue.colname],lambda=0)
    out.sites.adj <- cbind(sites.data, qvalue=q.sites)
    out.sites.adj.sel <- subset(out.sites.adj, subset=qvalue<cutoff)
  }
  else if(method =="BY" | method =="bonferroni") {
    adj.sites <- p.adjust(sites.data[,pvalue.colname], method=method)
    out.sites.adj <- cbind(out.sites, p.adj=adj.sites)
    out.sites.adj.sel <- subset(out.sites.adj, subset=p.adj<cutoff)
  }
  chrs <- unique(out.sites.adj.sel$chr)
  for (achr in chrs) {
    cat("*** ", as.character(achr), " ***", "\n")
    sites.sub <- subset(out.sites.adj.sel, subset=chr==achr)
    patts.sub <- unique(sites.sub$pattern)
    for (apatt in patts.sub) {
      cat(as.character(apatt),"\t")
      sites.sub.patt <- subset(sites.sub, subset=pattern==apatt)
      ctrs.sub <- unique(sites.sub.patt$contrast)
      for (actr in ctrs.sub) {
        cat(as.character(actr),"\t")
        ids <- subset(sites.sub.patt, subset=contrast==actr, select="ID", drop=T)
        cat(length(ids), "\t")
        out.reads.i <- subset(reads.data, subset=(chr==achr &
                                                  pattern == apatt &
                                                  contrast == actr &
                                                  ID %in% ids))
        out.reads.sel <- rbind(out.reads.sel, out.reads.i)
        cat(nrow(out.reads.i), "\n")
      }
    }
    
  }
  if(!is.null(sites.out.file)) {
    write.csv(out.sites.adj.sel, file=paste(out.path,
                                   paste(sites.out.file,"adj",method,cutoff,"csv",sep="."),
                                   sep="/"))
  }
  if(!is.null(reads.out.file)) {
    write.csv(out.reads.sel, file=paste(out.path,
                                 paste(reads.out.file,"adj",method,cutoff,"csv",sep="."),
                                 sep="/"))  
  }
  list(sites=out.sites.adj.sel,
       reads=out.reads.sel)
}
###.. Newer post-differential-test process
##--------------------------------------

select.wins.test.element <- function(x, method="fdr",
                                        cutoff=0.05,
                                        pvalue.col,
                                        effect.col,
                                        sign.names=c("Hypo","Hyper"))
{
  ##browser()
  pvalue.colname <- names(x)[pvalue.col]
  effect.colname <- names(x)[effect.col]
  
  sel <- subset(x, subset=p.adjust(x[,pvalue.colname], 'fdr') < 0.05)
  sign.idx <- (sel[,effect.colname] > 0) + 1
  output <- data.frame(sel, sign=sign.names[sign.idx])[,
                              c('start','end','ID',
                                "X.Intercept..Value",
                                effect.colname,
                                pvalue.colname,
                                "sign")]
  names(output) <- c('start','end','ID',"baseline","effect", "p.value", "sign")
  output
}

select.wins.test.report <- function(xl, method="fdr",
                                        cutoff=0.05,
                                        pvalue.col.list,
                                        effect.col.list,
                                        sign.names=c("Hypo","Hyper"))
{
  ## pvalue.col can be a vector: different elements in xl might take different p value cols.
  ##browser()
  n.list <- length(xl)
  output <- vector("list", length=n.list)

  ## if(length(pvalue.col.list) ==1) pvalue.col.list <- rep(pvalue.col.list, n.list)
  ## else stopifnot(length(pvalue.col.list) == n.list)

  ## if(length(effect.col) ==1) effect.col.list <- rep(effect.col.list, n.list)
  ## else stopifnot(length(effect.col.list) == n.list)
    
  for (i in seq(along = xl)) {
    pvalue.cols.i <- pvalue.col.list[[i]]
    effect.cols.i <- effect.col.list[[i]]
    output[[i]] <- vector("list", length(pvalue.cols.i))
    for (j in seq(along = pvalue.cols.i)) {
      
      pvalue.col.ij <- pvalue.cols.i[[j]]
      effect.col.ij <- effect.cols.i[[j]]
      xl.names.ij <- names(xl[[c(i,j)]])
      worker.ij <- function(k) {
        cat(i,"\t",j, "\t", k,"\n")
        ##pvalue.colname.j <- xl.names.ij[pvalue.col.j[k]]
        ##effect.colname.j <- xl.names.ij[effect.col.j[k]]
        ret <- select.wins.test.element(xl[[c(i,j)]],
                                        method = method,
                                        cutoff = cutoff,
                                        pvalue.col = pvalue.col.ij[k],
                                        effect.col = effect.col.ij[k],
                                        sign.names = sign.names)
        ret
      }
      output[[c(i,j)]] <- sapply(seq(along=pvalue.col.ij), worker.ij, simplify=F)
    }
  }
  output
}
  
    
  ## output <- lapply(wins.test.report,
  ##                  function(x) lapply(x,
  ##                                     function(y) select.wins.test.element(y,
  ##                                                                          method=method,
  ##                                                                          cutoff=cutoff,
  ##                                                                          pvalue.colname=pvalue.colname,
  ##                                                                          effect.colname=effect.colname,
  ##                                                                          sign.names=sign.names)))
##   output
## }

overlap.wins.test.element <- function(x1, x2, e, self=F) {
  poi.cols <- c('start','end')
  output <- any.overlap.sites(x1[,poi.cols],x2[,poi.cols],e=e,self=self)
  output
}

comb.hyperHypo <- function(wins.test.hyper,
                           wins.test.hypo,
                           v.varname,
                           p.varname) {
  wins.test.comb <- rbind(subset(wins.test.hyper, subset=wins.test.hyper[, v.varname] > 0),
                          subset(wins.test.hypo, subset=wins.test.hypo[,v.varname] < 0))
  wins.test.cutoff <- select.wins.test.element(wins.test.comb,
                                               pvalue.colname=p.varname,
                                               effect.colname=v.varname)
  wins.test.cutoff
}
 



###.. Clustering
##--------------------------------------
consolidate.lows <- function(dat, cols=2:3, limit=0, method=c("sum","max", "low")) {
  method <- match.arg(method)
  if(method == "sum") {
    idx <- rowSums(dat[,cols]) > limit
    
  }
  else if (method == "max") {
    idx <- apply(dat[, cols], 1, max) > limit
  }
  else if (method == "low") {
     idx <- apply(dat[, cols], 1, max) <= limit
  }
  else stop("Not implemented yet")
  res <- dat[idx,]
}

plot.multivariate <- function(dist.cor, which=c(1,2), clust.method="complete",...)
{ 
  ## 1 for cluster and 2 for mds
  if(1 %in% which) {
    clust <- hclust(dist.cor, method=clust.method)
    plot(as.dendrogram(clust),...)
    
  }
  if(2 %in% which) {
    mds <- cmdscale(dist.cor,k=4)
    plot(mds[,1], mds[,2], bty="n", xlab="Dimension 1", ylab="Dimension 2",...)
    text(mds[,1], mds[,2], labels=row.names(mds), pos=1, cex=0.7)
    plot(mds[,1], mds[,3], bty="n", xlab="Dimension 1", ylab="Dimension 3",...)
    text(mds[,1], mds[,3], labels=row.names(mds), pos=1, cex=0.7)
    plot(mds[,2], mds[,3], bty="n", xlab="Dimension 2", ylab="Dimension 3",...)
    text(mds[,2], mds[,3], labels=row.names(mds), pos=1, cex=0.7)
  }
}
###. process report (NEW)
## TODO: comment the above old process report code after the NEW works 
##____________________________________________________________________
###. Temerary file to contain functions used in processReport
##____________________________________________________________________

code.pattern <- function(patt.mx) {
  if(is.null(dim(patt.mx))) patt.mx <- as.matrix(patt.mx)
  cpatt <- apply(patt.mx, 2, paste, collapse="")
  gsub("NA","*", cpatt)
}

code.contrast <- function(contr.mx) {
  if(is.null(dim(contr.mx))) contr.mx <- as.matrix(contr.mx)
  colcontr <- apply(contr.mx, 2, paste, collapse="")
  paste(colcontr, collapse="_")
}

flattern.testReport <- function(test.report,
                                chr,
                                winsize,
                                colnames)
{
  ## test.report can either be wins.test or wins.test.report
  ##browser()
  ### Assume the same number of columns in the contrasts of test.report
  if(missing(colnames)) {
    ##n.ctr <- ncol(attr(test.report[[1]][[1]], "contrast")) (Not as reliable)
    n.ctr <- nrow(attr(test.report[[1]][[1]], "contrast"))-1
    seq.ctr <- seq(0,n.ctr)
    colnames=c("start","end","ID",
      paste("value",seq.ctr, sep="."),
      paste("p.value",seq.ctr, sep="."))
  }
  
  output <- NULL
  for(i in seq(along=test.report)) {
    for(j in seq(along=test.report[[i]])) {
      output.ij <- cbind(chr,
                         pattern=code.pattern(attr(test.report[[i]][[j]],'pattern')),
                         contrast=code.contrast(attr(test.report[[i]][[j]],'contrast')),
                         test.report[[i]][[j]]
                         )
      names(output.ij) <- c("chr","pattern","contrast",colnames)
      output <- rbind(output, output.ij)
    }
  }
  ###: CREATED LAST WIN PROBLEM output[,'end'] <- output[,'end']+winsize
  ##test.report.names <- names(test.report[[i]][[j]]) ## purposely
  ##names(output) <- c('chr','pattern','contrast',test.report.names)
  output
}

sel.sigsites <- function(flat.test.report,
                         pattern,
                         contrast,
                         method =c(
                           "holm",
                           "hochberg",
                           "hommel",
                           "bonferroni",
                           "BH",
                           "BY",
                           "fdr",
                           "none",
                           "qvalue"),
                         pvalue.var="p.value.1",
                         cutoff=0.05)
{
  ## See Table below from flat.test.report on columns pattern and contrast
  ## table(test.report.flat[,c('pattern','contrast')])
  ##          contrast
  ## pattern 011_0-11 001_-110 010_-101 010_001
  ##     011     5780        0        0       0
  ##     001        0      749        0       0
  ##     100     1783        0        0       0
  ##     110        0     4639        0       0
  ##     010        0        0      827       0
  ##     101        0        0     4112       0
  ##     111    30992    30992    30992   30992
  ##     **1        0        0        0   36847
  ##     *1*        0        0        0   36719
  ##     1**        0        0        0   36947

  ##browser()
  method <- match.arg(method)
  
  idx.patt <- flat.test.report$pattern %in% pattern
  idx.contr <- flat.test.report$contrast %in% contrast
  idx.start <- idx.patt & idx.contr

  idx.out <- NULL
  if(sum(idx.start)>0) {
    pvals <- subset(flat.test.report,
                    subset=idx.start,
                    select=pvalue.var,
                    drop=T)
    if(method != "qvalue") {
      pvals.adj <- p.adjust(pvals, method=method)
    } else if(method == "qvalue") {
      pvals.adj <- compute.qvalues(pvals,lambda=0)
    } else stop("Not supported multiple testing method provided!")
    
    idx.out <- which(idx.start)[pvals.adj < cutoff]
  }
  idx.out
}

compute.qvalues <- function(pvalues,...) {
  qvalue::qvalue(pvalues, gui=F,...)$qvalues
}

select.sig.sites <- function(flat.test.report,
                             method.p.sel="qvalue",
                             cutoff.p.sel=0.05,
                             cutoff.v.sel=NA,
                             winsize=50,
                             cutoff.by.winsize=T,
                             cutoff.by.winsize.pct=0.9,
                             pvalue.var="p.value.1",
                             value.var="value.1",
                             split.vars=c("pattern", "contrast"),
                             cols.remove=c("space","seq", "ID")) {
  
  method.p.sel <- match.arg(method.p.sel)
  if(method.p.sel == "qvalue") {
    .rd.full <- trans2rangedData(dpt.report.postprocess(flat.test.report,
                                                        cutoff.by.winsize=cutoff.by.winsize,
                                                        winsize=winsize,
                                                        cutoff.by.winsize.pct=cutoff.by.winsize.pct,
                                                        pvalue.var="p.value.1",
                                                        ##value.var="value.1",
                                                        split.vars=split.vars,
                                                        cols.remove=cols.remove),
                                 ext.end=NULL)
    
    .rd.sel <- cutoff.sites(.rd.full,
                            qvalue.cutoff=cutoff.p.sel,
                            value.cutoff=cutoff.v.sel,
                            value.var=value.var,
                            qvalue.var="qvalue")
    .df.sel <- as.data.frame(.rd.sel)
    .df.sel <- .df.sel[,!(names(.df.sel)%in%c("space"))]
    
  } else {stop("Cutoff method not implemented!")}
  
  attr(.df.sel, 'testing') <- list(method=method.p.sel, cutoff=cutoff.p.sel)
  .df.sel
}

replace.names <- function(names, target, replace) {  
  names[match(target, names)] <- replace
  names
}
           
cutoff.sites <- function(q.dt,
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

dpt.report.postprocess <- function(flat.test.report,
                                   cutoff.by.winsize=T,
                                   winsize=50,
                                   cutoff.by.winsize.pct=0.9,
                                   pvalue.var="p.value.1",
                                   split.vars=c("pattern", "contrast"),
                                   cols.remove=c("space","seq", "ID")) {
  ##browser()
  .dt <- trans2rangedData(remove.dupsite(flat.test.report), ext.end=winsize-1)
  .dt$qvalue <- dpt.qvalues(.dt,
                            cutoff.by.winsize=cutoff.by.winsize,
                            winsize=winsize,
                            cutoff.by.winsize.pct=cutoff.by.winsize.pct,
                            pvalue.var=pvalue.var,
                            split.vars=split.vars)
  .dt <- remove.dupcontr(.dt, pvalue.var=pvalue.var)
  .dt[,!(names(.dt)%in%cols.remove)]
}

trans2rangedData <- function(df.in, ext.end=49) {
  ## df.in needs to have start, end, and chr as columns
  ## end in df.in is the start position of the last window
  ## ext.end set to NULL for properlly calculated ends.
  if(!is.null(ext.end)) df.in$end <- df.in$end+ext.end
  df.in <- cbind(df.in, space = df.in$chr)
  as(df.in, "RangedData")
}

dpt.qvalues <- function(flat.test.report,
                        cutoff.by.winsize=T,
                        winsize=50,
                        cutoff.by.winsize.pct=0.9,
                        pvalue.var="p.value.1",
                        split.vars=c("pattern", "contrast")) {
  
### handle NA pvalues
### flat.test.report is a RangedData class
  qvalues.init <- rep(NA, nrow(flat.test.report))
  pval.na.idx <- is.na(flat.test.report[[pvalue.var]])  

### split by pattern, contrast, and possibly winsize
  cat("# of sites:", nrow(flat.test.report), "\n")
  cat("# of NA pvalues:", sum(pval.na.idx), "\n")
  
  in.sub <- flat.test.report[!pval.na.idx,]
  pvals <- in.sub[[pvalue.var]]
  in.width <- width(in.sub)
  
  fac.list <- list(patt=in.sub$pattern, ctr=in.sub$contrast)

  pvals.split <- split(data.frame(p=pvals, w=in.width), fac.list)
  qvals.split <- lapply(pvals.split,
                        function(x, cut.pct) {
                         
                          if(nrow(x)==0) {
                            return(numeric())
                          } else {
                            w <- x$w
                            p <- x$p
                            w.fac.max <- quantile(w, cut.pct)
                            w[w>w.fac.max] <- w.fac.max+1
                            w.fac <- factor(w)
                            p.split <- split(p, w.fac)
                            unsplit(lapply(p.split, compute.qvalues, lambda=0), w.fac)
                          }
                        },
                        cut.pct=cutoff.by.winsize.pct)
  
### unsplit the calculated qvalues
  qvalues <- unsplit(qvals.split, fac.list)
  qvalues.init[!pval.na.idx] <- qvalues
  qvalues.init
}

remove.dupsite <- function(flat.test.report) {
### flat.test.report is a data.frame
  unique(flat.test.report)
}

remove.dupcontr <- function(flat.test.report.afterqvalue,
                            pvalue.var="p.value.1") {
  ##browser()
### flat.test.report.afterqvalue is a RangedData class or data.frame
  .df <- as.data.frame(flat.test.report.afterqvalue[!is.na(flat.test.report.afterqvalue[[pvalue.var]]),])
  .df$seq <- seq(nrow(.df))
  .df$id <- with(.df, paste(chr, pattern, ID, sep="_"))
  seq.select <-  which.min.by(.df$id, .df[[pvalue.var]])
  .df[seq.select,]
  
}

which.min.by <- function(by, score) {
  .id <- seq(length(score))
  score.order <- order(score, na.last=T,
                       decreasing=F)
  .select <- vaggregate(.id[score.order],
                        factor(by[score.order]),
                        function(x) x[1])
  .select
}


which.max.by <- function(by, score) {
  
  .id <- seq(length(score))
  score.order <- order(score, na.last=T,
                       decreasing=T)
  .select <- vaggregate(.id[score.order],
                        factor(by[score.order]),
                        function(x) x[1])
  .select
}


select.sig.sites.2 <- function(flat.test.report, contrast.objs, method.sel="fdr", cutoff.p.sel=0.05,
                               winsize=50, cutoff.by.winsize=T, ul.by.winsize=30)
{

  ### remove NA pvalues
  pval.na.idx <- which(is.na(flat.test.report$p.value.1))
  if(length(pval.na.idx) > 0) {
    cat("# of missing p values: ", length(pval.na.idx),"\n")
    flat.test.report <- flat.test.report[-pval.na.idx,]
  }
  ### process selection
  
  out <- NULL
  if(cutoff.by.winsize) {
    wsizes <- as.integer(names(with(flat.test.report, table(end-start))))
    wsize.ul <- ul.by.winsize * winsize
    singles <- wsizes[wsizes<=wsize.ul]
    for(asize in singles) {
      a.report <- subset(flat.test.report, subset=(end-start)==asize)
      idx.sel <- NULL
      for(i in seq(along=contrast.objs)) {
        idx.sel <- union(idx.sel, sel.sigsites(a.report,
                                               pattern = contrast.objs[[i]]$pattern,
                                               contrast = contrast.objs[[i]]$contrast,
                                               method = method.sel,
                                               pvalue.var = contrast.objs[[i]]$pvalue.var,
                                               cutoff = cutoff.p.sel)
                         )
      }
      out <- rbind(out, a.report[idx.sel,])
    }
    a.report <- subset(flat.test.report, subset=(end-start)>wsize.ul)
    idx.sel <- NULL
    for(i in seq(along=contrast.objs)) {
      idx.sel <- union(idx.sel, sel.sigsites(a.report,
                                             pattern = contrast.objs[[i]]$pattern,
                                             contrast = contrast.objs[[i]]$contrast,
                                             method = method.sel,
                                             pvalue.var = contrast.objs[[i]]$pvalue.var,
                                             cutoff = cutoff.p.sel)
                       )
    }
    out <- rbind(out, a.report[idx.sel,])
  } else {
    idx.sel <- NULL
    for(i in seq(along=contrast.objs)) {
      idx.sel <- union(idx.sel, sel.sigsites(flat.test.report,
                                             pattern = contrast.objs[[i]]$pattern,
                                             contrast = contrast.objs[[i]]$contrast,
                                             method = method.sel,
                                             pvalue.var = contrast.objs[[i]]$pvalue.var,
                                             cutoff = cutoff.p.sel)
                       )
    }
    out <- test.report.flat[idx.sel,]
  }
  attr(out, 'contrasts') <- contrast.objs
  attr(out, 'testing') <- list(method=method.sel, cutoff=cutoff.p.sel)
  out
}
                                           

  
append.effectsign <- function(test.report,
                             col="value.1",
                             colname="MethylType",
                             signvalues=c("Hyper", "Hypo")
                             )
{
  signs <- test.report[,"value.1"] >0
  signvalues <- rev(signvalues)[signs+1]
  out <- cbind(test.report, signvalues)
  names(out) <- c(names(test.report), colname)
  out
}

get.difftest.flatReport <- function(path,
                                    chrs,
                                    winsize,
                                    pos=-1
                                    )
{
  ##browser()
  dpt.options <- get("dpt.options", pos=pos)
  ws <- 'diffTest'
  if(missing(path)) path <- get.ws.path(ws, pos=pos)
  if(missing(chrs)) chrs <- get.ws.chrs(ws, pos=pos)
  report.flat <- NULL
  for (chr in chrs) {
    cat("loading",dpt.options[[c("output","diffTest")]],"data:",chr)
    load.diffTestData(path, chr, pos=pos)
    cat("\t done\n")
    report.flat <- rbind(report.flat,
                         flattern.testReport(wins.test.report,
                                             chr,
                                             winsize))
  }
  report.flat
}

###. Workspace (ws) utilities
##____________________________________________________________________

## Assumption: session is running within the workspace directory.

get.wsoption.path <- function(which=c("root",
                                "mixPoisson",
                                "pattRecog",
                                "diffTest",
                                "preprocess",
                                "joinSample",
                                "log",
                                "postprocess",
                                "report",
                                "type",
                                "analysis",
                                "src",
                                "tmp"),
                              pos=-1
                              ## those which options are mapped from all possible places
                              ## that might get queried to workspace
                              )
{
  which <- match.arg(which)
  dpt.options <- get("dpt.options", pos=pos)
  dpt.ws.options <- dpt.options[['ws']]
  if(which == "root") dpt.ws.options$ws.root
  else if(which %in% c("type","analysis")) dpt.ws.options$analysis.path
  else if(which == "mixPoisson") dpt.ws.options$mixpos.path
  else if(which == "pattRecog") dpt.ws.options$pattrecog.path
  else if(which == "diffTest") dpt.ws.options$difftest.path
  else if(which %in% c("preprocess","joinSample")) dpt.ws.options$preprocess.path
  else if(which == "log") dpt.ws.options$log.path
  else if(which %in% c("postprocess","report")) dpt.ws.options$postprocess.path
  else if(which == "src") dpt.ws.options$src.path
  else if(which =="tmp") dpt.ws.options$tmp.path
  else stop("Not a valid option: ", which)
}

get.ws.path <- function(which=c("root",
                          "mixPoisson",
                          "pattRecog",
                          "diffTest",
                          "joinSample",
                          "log",
                          "report",
                          "analysis",
                          "src",
                          "tmp"),
                        pos=-1
                        ## those which options are very specificly matched to tasks
                        )
{
  which <- match.arg(which)
  dpt.options <- get("dpt.options", pos=pos)
  if(file.exists(get.wsoption.path("root", pos=pos))) {
    ws.root <- get.wsoption.path("root", pos=pos)
  } else if(file.exists(paste("..",get.wsoption.path("root", pos=pos), sep="/"))) {
    ws.root <- paste("..",get.wsoption.path("root",pos=pos), sep="/")
  } else {
    ws.root <- file.path(dpt.options[[c("ws","out.path")]],
                         get.wsoption.path("root", pos=pos))
  }
  if(which == "root") ws.root
  else file.path(ws.root, get.wsoption.path(which, pos=pos))
}


create.ws.dir <- function(which=c("root",
                                "mixPoisson",
                                "pattRecog",
                                "diffTest",
                                "preprocess",
                                "log",
                                "report",
                                "src",
                                "tmp"),
                          pos=-1
                          )
{
  which <- match.arg(which)
  .dir <- get.wsoption.path(which, pos=pos)
  if(!file_test("-d", .dir)) {
    .create.yes <- dir.create(.dir)
    if(!.create.yes) stop(paste("Can not create",
                                .dir,
                                "for dpt analyses!",
                                sep=" "))
  } else {
    cat(.dir, "exists!","\n")
  }
}

init.ws.dirs <- function(what=c("log","src","tmp"), pos=-1) {
  create.ws.dir("root")
  ws.root.path <- get.ws.path('root', pos=pos)
  current.dir <- getwd()
  setwd(ws.root.path)
  for (ipath in what) {
    create.ws.dir(ipath, pos=pos)
  }
  setwd(current.dir)
}

save.ws.image <- function(ws.root.path, pos=-1) {
  if(missing(ws.root.path)) ws.root.path <- get.ws.path('root',pos=pos)
  current.dir <- getwd()
  setwd(ws.root.path)
  save.image()
  setwd(current.dir)
}

load.ws.image <- function(ws.root.path, pos=-1) {
  if(missing(ws.root.path)) ws.root.path <- get.ws.path('root', pos=pos)
  current.dir <- getwd()
  setwd(ws.root.path)
  load(".RData"); .First()
  setwd(current.dir)
}

get.wsoption.para <- function(which=c("winsize"), pos=-1) {
  which <- match.arg(which)
  if(which == "winsize") get.winsize()
  else stop("Not a valid DiffSeq parameter: ", which)
}

get.n.cpus <- function(n.cpus) {
  ncpus <- as.numeric(gsub("cpus","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) ==0) ncpus <- as.numeric(gsub("ncpus","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) ==0) ncpus <- as.numeric(gsub("ncores","",gsub('-','',grep('[:digit:]*cpus', commandArgs(), v=T))))
  if(length(ncpus) == 0) ncpus <- 1
  if(missing(n.cpus)) ncpus else n.cpus
}

get.commandArg <- function(key, mode=c("numeric","character")) {
  mode <- match.arg(mode)
  if(mode == "numeric") argInput <- as.numeric(gsub(key,"", gsub('-','',grep(key, commandArgs(), v=T))))
  else argInput <- gsub(key,"", gsub('-','',grep(key, commandArgs(), v=T)))
  if(length(argInput)==0) stop("Please provide a command argument as '--cutoff-0.01'!")
  argInput
}
  

get.winsize.input <- function(winsize) {
  wsize <- as.numeric(gsub("winsize","",
                           gsub('-','',grep('winsize[:digit:]*', commandArgs(), v=T))))
  if(length(wsize) == 0) wsize <- 100
  if(missing(winsize)) wsize else winsize
}

get.winsize <- function(path, pos=-1) {
  if(missing(path)) winsize.path <- get.ws.path("joinSample", pos=pos)
  load(paste(winsize.path, "winsize.RData",sep="/"))
  winsize
}


get.str.chr <- function(chrStr) {
  chr <- gsub('-','',grep('chr', commandArgs(), v=T))
  if(missing(chrStr)) {
    if(length(chr) == 0) stop("chr string is missing.")
    else chr
  } else chrStr
}


get.ws <- function(ws=c("mixPoisson","pattRecog","diffTest","joinSample","log","report","analysis"),
                   what,
                   chr,
                   ws.path,
                   pos=parent.frame(),
                   pos.out=parent.frame())
{
  ##output.list <- NULL
  ws <- match.arg(ws)  
  ##if(missing(pos.out)) pos.out <- parent.frame()
  if(ws == "mixPoisson") {
    if(missing(chr)) chr <- get('chr', pos=pos)
    if(missing(ws.path)) ws.path <- get.ws.path("mixPoisson", pos=pos)
    load.mixPoiData(ws.path, chr, pos=pos)
    for (eachwhat in what) {
      ##output.list <- list(output.list, get(eachwhat))
      assign(eachwhat, get(eachwhat), pos=pos.out) 
    }
  }
  else if(ws == "pattRecog") {
    if(missing(chr)) chr <- get('chr', pos=pos)
    if(missing(ws.path)) ws.path <- get.ws.path("pattRecog", pos=pos)
    load.pattRecogData(ws.path, chr, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
  }
  else if(ws == "diffTest") {
   if(missing(chr)) chr <- get('chr', pos=pos)
   if(missing(ws.path)) ws.path <- get.ws.path("diffTest", pos=pos)
   load.diffTestData(ws.path, chr, pos=pos)
   for (eachwhat in what) {
     assign(eachwhat, get(eachwhat), pos=pos.out)
   }
 }
  else if(ws == "report") {
    if(missing(ws.path)) ws.path <- get.ws.path("report", pos=pos)
    load.reportData(ws.path, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
    
  }
  else if(ws == "analysis") {
    if(missing(ws.path)) ws.path <- get.ws.path("report", pos=pos)
    load.analysisData(ws.path, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
    
  }
  ##names(output.list) <- what
  ##output.list
}

get.reads.fromJoinedSamples <- function(inputPath, chr, threshold=1, pos=-1,
                                        pos.out=parent.frame(), aug.fun=max,
                                        sample.select=NULL) {
  dpt.options <- get("dpt.options", pos=pos)
  reads.datafile <- paste(inputPath,
                          paste(dpt.options[[c("output","joinSample")]],chr,"RData", sep="."),
                          sep="/")
  lf <- load(reads.datafile)
  
  s.cols <- grep("^s[0-9]*",names(get(lf)))
  if(!is.null(sample.select)) s.cols <- s.cols[sample.select]
  n.sample <- length(s.cols)
  reads <- get(lf)[apply(get(lf)[,s.cols,drop=F], 1, aug.fun) > threshold,]
  reads.poi <- reads$poi
  reads.names <- names(reads)[s.cols]
  assign("s.cols", s.cols, pos = pos.out)
  assign("reads", reads, pos=pos.out)
  assign("reads.poi", reads.poi, pos=pos.out)
  assign("reads.names", reads.names, pos=pos.out)
}

get.mappedCountTab <- function(inPath=NULL, pos=-1) {
  if(is.null(inPath)) {
    inPath <- get.wsoption.path("joinSample", pos=pos)
  }
  cntab.data.file <- file.path(inPath, "countdistr.RData")
  load(cntab.data.file)
  ##assign("counts.tab", count.tbs, pos=parent.frame(n=1))
  lapply(count.tbs, function(x) {x$count <- as.numeric(as.character(x$count));x})
  ##count.tbs
}
    
compute.pmeans <- function(res, reads.poi, reads.names, count.table, pos.out=parent.frame()) {

  ## if(!missing(ws.mixPois.path)) load.mixPoiData(ws.mixPois.path)
  ## res <- res[!unlist(lapply(res, is.null))]
  
  pmean.res <- as.data.frame(lapply(res, get.postMean, log=F))
  row.names(pmean.res) <- reads.poi
  names(pmean.res) <- reads.names

  pmean.var <- as.data.frame(lapply(res, get.postmean.var))
  row.names(pmean.var) <- reads.poi
  names(pmean.var) <- reads.names

  pmeans <- list(##lambda2 = lambda2.res,
                 mean = pmean.res,
                 var = pmean.var)
  
  ##pmeans.norm <- normalize.sig(pmeans$mean, pmeans$lambda2)
  ##pmeans.norm <- normalize.pmeans(pmeans$mean, countab)
  pmeans.norm <- as.data.frame(normalize.quantiles(as.matrix(pmean.res)))
  row.names(pmeans.norm) <- reads.poi
  names(pmeans.norm) <- reads.names

  ### bg corrected
  pmeans.bgcorrect.norm <- as.data.frame(normalize.quantiles(as.matrix(normalize.pmeans(pmean.res,
                                                                                        count.table,
                                                                                        bg.correct=T,
                                                                                        out.int=F))))
  row.names(pmeans.bgcorrect.norm) <- reads.poi
  names(pmeans.bgcorrect.norm) <- reads.names

  ### bg normalized but not corrected
  pmeans.bgnorm.norm <- as.data.frame(normalize.quantiles(as.matrix(normalize.pmeans(pmean.res,
                                                                                        count.table,
                                                                                        bg.correct=F,
                                                                                        out.int=F))))
  row.names(pmeans.bgnorm.norm) <- reads.poi
  names(pmeans.bgnorm.norm) <- reads.names

  assign("pmeans", pmeans, pos = pos.out)
  assign("pmeans.norm", pmeans.norm, pos = pos.out)
  assign("pmeans.bgnorm.norm", pmeans.bgnorm.norm, pos = pos.out)
  assign("pmeans.bgcorrect.norm", pmeans.bgcorrect.norm, pos = pos.out)
  
}


## get.sites <- function(ws.pattRecog.path,
##                       chr,
##                       which=c("sites",
##                         "sites.j",
##                         "sites.j.ex",
##                         "sites.js",
##                         "sites.js.ex"),
##                       pos=-1,
##                       pos.out=parent.frame())
## {
##   ##?? Not working, use get.ws instead.
##   if(missing(ws.pattRecog.path)) ws.pattRecog.path <- get.ws.path("pattRecog", pos=pos)
##   ##chr.env <- parent.env(environment())
##   if(missing(chr)) chr <- get("chr", pos=pos, inherits = TRUE)
##   load(file.path(ws.pattRecog.path, paste("pattRecog_wins",chr,"RData", sep=".")))
##   for(eachwhich in which) {
##     assign(eachwhich, get(eachwhich), pos=pos)
##   }
## }

## save.ws.diffTest <- function(outPath, what, chr, file.header,pf.n=1) {
##   if(missing(outPath)) outPath <- get.ws.path("diffTest")
##   ## chr.env <- parent.env(environment())
##   ## if(missing(chr)) chr <- get("chr", envir=chr.env, inherits = TRUE)
##   save(list=what,
##        file =file.path(outPath, paste(file.header, chr, "RData", sep=".")),
##        envir=parent.frame(n=pf.n))
## }


## save.ws.mixPoisson <- function(outPath, what, chr, file.header, pf.n=1) {
##   if(missing(outPath)) outPath <- get.ws.path("mixPoisson")
##   save(list=what,
##        file=file.path(outPath, paste(file.header,chr,"RData", sep=".")),
##        envir=parent.frame(n=pf.n)
##        )
## }

## save.ws.pattRecog <- function(outPath, what, chr, pf.n=1) {
##   if(missing(outPath)) outPath <- get.ws.path("pattRecog")
##   save(list=what,
##        file=file.path(outPath, paste("pattRecog",chr,"RData", sep=".")),
##        envir=parent.frame(n=pf.n)
##        )
## }

## save.ws.report <- function(outPath, what, pf.n=1) {
##   if(missing(outPath)) outPath <- get.ws.path("report")
##   save(list=what,
##        file=file.path(outPath, paste("report","RData", sep=".")),
##        envir=parent.frame(n=pf.n)
##        )
## }

## save.ws.joinSamples <- function(outPath, what, chr, pf.n=1) {
##   if(missing(outPath)) outPath <- get.ws.path("joinSamples")
##   save(list=what,
##        file=file.path(outPath, paste("joinedSamples",chr,"RData", sep=".")),
##        envir=parent.frame(n=pf.n)
##        )
## }
  
save.ws <- function(ws=names(dpt.options[['output']]),
                    what,
                    outpath,
                    chr=NULL,
                    pos=parent.frame()) {
  ws <- match.arg(ws)
  dpt.options <- get("dpt.options", pos=pos)
  if(missing(outpath)) outpath <- get.ws.path(ws, pos=pos)
  if(is.null(chr)) {
    outfile <- file.path(outpath, paste(dpt.options[[c('output',ws)]],"RData", sep="."))
  } else {
    outfile <- file.path(outpath, paste(dpt.options[[c('output',ws)]],chr,"RData", sep="."))
  }
  ##save.envir <- if(pos==-1) globalenv() else as.environment(pos)
  save(list=what,
       file=outfile,
       envir=as.environment(pos))
  cat(what,"from",ws,"\t","saved!\n")
}
    
get.str.smpsheet <- function(smpsheetStr) {
  smpsheet <- gsub('-','',grep('samplesheet', commandArgs(), v=T))
  if(missing(smpsheetStr)) smpsheet else smpsheetStr
}

load.mappedData <- function(smpStr, from=c("R","file"), pos=-1, pos.out=parent.frame()) {
  ##browser()
  from <- match.arg(from)
  if(from == "R") infiles <- smpStr
  else if(from == "file") infiles <- scan(get.str.smpsheet(smpStr), what="list", comment.char="#")
  files <- paste("..",infiles,sep="/")
  dataload <- list()

  for (i in seq(along=files)) {
    dataload <- c(dataload, list(read.table(files[i])))
  }
  names(dataload) <- paste("s", seq(along=dataload), sep="")
  n.reads <- unlist(lapply(dataload, function(x) sum(x$V4)))
  save(n.reads, file=file.path(get.wsoption.path("joinSample", pos=pos),"nreads.RData"))

  count.tbs <- sapply(dataload, function(x) {idf <- as.data.frame(table(x$V4));
                                             names(idf) <- c("count","freq");
                                             idf},
                      simplify=F)
  save(count.tbs, file=file.path(get.wsoption.path("joinSample", pos=pos),"countdistr.RData"))
  
  assign("dataload", dataload, pos=pos.out)
  assign("n.reads", n.reads, pos=pos.out)
}

load.mixPoiData <- function(from, chr, pos=-1, pos.out=parent.frame()) {
  if(missing(chr)) chr <- get('chr', pos=pos)
  if(missing(from)) from <- get.ws.path('mixPoisson', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from,paste(dpt.options[[c('output','mixPoisson')]],chr,"RData", sep=".")),
       envir=as.environment(pos.out))
}

load.pattRecogData <- function(from, chr, pos=-1, pos.out=parent.frame()) {
  if(missing(chr)) chr <- get('chr', pos=pos)
  if(missing(from)) from <- get.ws.path('pattRecog', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','pattRecog')]],chr,"RData", sep=".")),
       envir=as.environment(pos.out))
}

load.diffTestData <- function(from, chr, pos=-1, pos.out=parent.frame()) {
  if(missing(chr)) chr <- get('chr', pos=pos)
  if(missing(from)) from <- get.ws.path('diffTest', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','diffTest')]],chr,"RData", sep=".")),
       envir=as.environment(pos.out))
}

load.reportData <- function(from, pos=-1, pos.out=parent.frame()) {
  if(missing(from)) from <- get.ws.path('report', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','report')]],"RData", sep=".")),
       envir=as.environment(pos.out))
}

load.analysisData <- function(from, pos=-1,pos.out=parent.frame()) {
  if(missing(from)) from <- get.ws.path('report', pos=pos)
  dpt.options <- get("dpt.options", pos=pos)
  load(file.path(from, paste(dpt.options[[c('output','analysis')]],"RData", sep=".")),
       envir=as.environment(pos.out))
}

get.ws.chrs <- function(which=c("mixPoisson",
                          "pattRecog",
                          "diffTest",
                          "joinSample",
                          "full"
                          ),
                        pos=-1
                        )
{
  which <- match.arg(which)
  dpt.options <- get("dpt.options", pos)
  if(which != "full") {
    files <- list.files(get.ws.path(which, pos=pos))
    if(which == "mixPoisson") {
      filestrs <- grep(paste(dpt.options[[c('output','mixPoisson')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "pattRecog") {
      filestrs <- grep(paste(dpt.options[[c('output','pattRecog')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "diffTest") {
      filestrs <- grep(paste(dpt.options[[c('output','diffTest')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which == "joinSample") {
      filestrs <- grep(paste(dpt.options[[c('output','joinSample')]],".chr[0-9A-Za-z]*.RData",sep=""),
                       files,
                       v=T)
    } else if(which =="full") {
      filestrs <- NULL
    } else stop("Not a valid option: ", which)
    out <- sapply(filestrs,
                  function(x) {
                    this.str <- unlist((strsplit(x, split="\\.")))
                    this.str[length(this.str)-1]
                  }
                  )
  }
  if(which == "full") out <- c(
       "chr10","chr11","chr12","chr13","chr14",
       "chr15","chr16","chr17","chr18","chr19",
       "chr1","chr20","chr21","chr22","chr2",
       "chr3","chr4","chr5","chr6","chr7","chr8",
       "chr9","chrX","chrY"
       )
  as.vector(out)
}

get.sharedSource <- function(url, pos=-1) {
  src <- getURL(url)
  ##ws.env <- if(pos==-1) globalenv() else as.environment(pos)
  eval(parse(text = src), as.environment(pos))
}

get.sharedSrc <- function(str, ver="current", pos=-1) {
  pathstr <- paste("http://dl.dropbox.com/u/10421230/dpt", ver, sep="/")
  urlStr <- paste(pathstr, str, sep="")
  get.sharedSource(urlStr, pos=pos)
}

## get.src <- function(string, base, dpt.option){
##   if(missing(dpt.option)) {
##     dpt.option <- list(ws.root = "DPT_ws",
##                        mixpos.path = "MixPoisson",
##                        pattrecog.path = "PattRecog",
##                        difftest.path = "DiffTest",
##                        joinsample.path = "JoinSamples",
##                        log.path = "logs",
##                        report.path = "Reports",
##                        src.path = "src")
##   }
##   attach(dpt.option)
##   if(!missing(base)) {
##     source(file.path(base, string))
##   } else {
##     if(file.exists(file.path(src.path, string))){
##       source(file.path(src.path, string))
##     } else if(file.exists(file.path(ws.root,src.path, string))) {
##       source(file.path(ws.root, src.path, string))
##     } else {stop("Can't allocate the src file")}
##   }
## }

match.patt.pattmx <- function(patt, pattmx) {
  match(patt, paste("P", apply(pattmx,2,paste, collapse=""), sep=""))
}

load.init.ws <- function(path, pos=-1) {
  if(missing(path)) {
    load(file.path(get.ws.path("root", pos=pos),".RData"))
  }
  else load(file.path(path, ".RData"))
}

###. Parallel DPT
##____________________________________________________________________
if.parallel <- function(ncpus) {
  if(length(ncpus)<1 || ncpus < 2) para.nostart=T else para.nostart=F
  para.nostart
}

start.para <- function(ncore, initExpr,...) {
  sfInit(parallel=TRUE, cpus=ncore, nostart=if.parallel(ncore),...)
  sfExport("initExpr")
  .dump <- sfClusterEval(eval(initExpr))
}

stop.para <- function(...) {
  sfStop(...)
}

paraApply <- function(x, margin, fun, ...,mode=real64(),mode.size=8, tmp.dir=".tmp",ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Memory efficient with mmap
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  if(missing(mode.size)) mode.size <- sizeof(mode)
  tmp <- tempfile(tmpdir=tmp.dir)
  xdim <- dim(x)
  .cut <- cut(seq(xdim[margin]), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  
  if(margin==1) {
    writeBin(as.vector(t(x)), tmp)
    frac.len <- ncol(x)
  } else {
    writeBin(as.vector(x), tmp)
    frac.len <- nrow(x)
  }
  
  fun.text <- deparse(substitute(fun))
  ##in.mmap <- mmap(tmp, mode=mode, )

  .worker <- function(y,
                      margin,
                      .file.mmap,
                      frac.len,
                      mode,
                      mode.size,
                      xdim,
                      ncore,
                      fun.text,
                      .cut,
                      ...) {
    ##browser()
    fun <- eval(parse(text=fun.text))
    .len <- sum(.cut==levels(.cut)[y])
    if(y>1) {
      .offset <- sum(.cut%in%levels(.cut)[1:(y-1)])
    } else {
      .offset <- 0
    }
    .data.mmap <- mmap(.file.mmap, mode=mode, prot=mmapFlags("PROT_READ"))
    ##cat(class(length(.data.mmap)))
    ##cat("which=",which)
    ovec <- tapply(.data.mmap[seq(.len*frac.len)+.offset*frac.len],
                   list(rep(seq(.len), each=frac.len)), fun,...)
    dimnames(ovec) <- NULL
    munmap(.data.mmap)
    ##cat(str(ovec))
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    margin=margin,
                    .file.mmap=tmp,
                    frac.len=frac.len,
                    mode=mode,
                    mode.size=mode.size,
                    xdim=xdim,
                    ncore=ncore,
                    fun.text=fun.text,
                    .cut=.cut,
                    simplify=F,
                    ...)
  ##munmap(in.mmap)
  unlink(tmp)
  unsplit(out.l, .cut)
}

paraApply2 <- function(x, margin, fun, ..., tmp.dir=".tmp",ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Higher throughput with shared data file
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  tmp <- tempfile(tmpdir=tmp.dir)
  xdim <- dim(x)
  .cut <- cut(seq(xdim[margin]), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  fun.text <- deparse(substitute(fun))
  save(x, margin, xdim, ncore, fun.text, .cut, file=tmp)
 
  .worker <- function(y,load.file,...) {
    ##browser()
    load(load.file)
    fun <- eval(parse(text=fun.text))
    .len <- sum(.cut==levels(.cut)[y])
    if(y>1) {
      .offset <- sum(.cut%in%levels(.cut)[1:(y-1)])
    } else {
      .offset <- 0
    }
    ovec <- apply(x[seq(.len)+.offset,],margin, fun,...)
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    load.file=tmp,
                    simplify=F,
                    ...)
  unlink(tmp)
  unsplit(out.l, .cut)
}

paraLapply <- function(x, fun, ..., tmp.dir=".tmp", ncore=NULL) {
  ## Parallel version apply with ncore blocks on the input matrix x
  ## Higher throughput with shared data file
  ##browser()
  if(is.null(ncore)) ncore <- sfCpus()
  if(is.null(tmp.dir)) tmp.dir <- get.ws.path("tmp")
  tmp <- tempfile(tmpdir=tmp.dir)
  xn <- length(x)
  .cut <- cut(seq(xn), breaks=ncore)
  levels(.cut) <- seq(along=levels(.cut))
  fun.text <- deparse(substitute(fun))
  save(x, xn, ncore, fun.text, .cut, file=tmp)
 
  .worker <- function(y,load.file,...) {
    ##browser()
    load(load.file)
    fun <- eval(parse(text=fun.text))
    .len <- sum(.cut==levels(.cut)[y])
    if(y>1) {
      .offset <- sum(.cut%in%levels(.cut)[1:(y-1)])
    } else {
      .offset <- 0
    }
    ovec <- lapply(x[seq(.len)+.offset,], fun,...)
    ovec
  }
  out.l <- sfSapply(seq(ncore),
                    .worker,
                    load.file=tmp,
                    ...)
  unlink(tmp)
  out.l
}

###... DPT process manager
batch.jobs <- function(caller, nc, nb, chrs=NULL, by.chr=T, retry.n=3, waiting=90) {
  .failed <- NULL
  ncore <- nc
  if(!by.chr) {
    .failed <- eval(caller)
  } else {
    .continue <- T
    .ntry <- 1
    while(.continue) {
      if(.ntry>1) cat("try = ", .ntry, "\n")
      .jobs <- (seq(along=chrs)-1)%/%nb
      .failed.try <- NULL
      for (i in unique(.jobs)) {
        .chrs <- chrs[.jobs==i]
        for (chr in .chrs){
          mcparallel(caller)
          dpt.syshold(waiting)
        }
        .fl <- c(unlist(collect()))
        .failed.try <- c(.failed.try, .fl)
      }
      chrs <- .failed.try
      .ntry <- .ntry+1
      .continue <- (.ntry<=retry.n)&(length(.failed.try)>0)
    }
    .failed <- .failed.try
  }
  return(.failed)
}    

dpt.syshold <- function(tm, verbose=F) {
  t1 <- proc.time()
  t.lag <- 0
  while(t.lag<tm) {
    system.time(for(i in 1:10000) x <- mean(rt(1000, df=4)))
    t2 <- proc.time()
    t.lag <- (t2-t1)[3]
  }
  if(verbose) cat("system waited for",t.lag,"sec\n")
}
###. Mis.
##____________________________________________________________________

siteUID <- function(tst) {
  paste(tst$chr,
        tst$pattern,
        tst$ID, sep="-")
}

clean.dupcontr.sites <- function(sites,
                                 p.value.col,
                                 n.para=10)
{
  ## this seems to be the best I can do in R
  ##browser()
  siteIDs <- siteUID(sites)
  cat("# of site IDs:", length(siteIDs),"\n")
  cat("# of unique site IDs:", length(unique(siteIDs)),"\n")
  dupIdx <- duplicated(siteIDs)
  #dupIDs <- unique(siteIDs[table(siteIDs)>1])
  dupIDs <- unique(siteIDs[dupIdx])
  cat("# of unique site IDs with duplicates:", length(dupIDs), "\n")
  idx.dup <- siteIDs%in%dupIDs
  cat("# of sites with duplicates", sum(idx.dup), "\n")
  out.0 <- data.frame(.id=siteIDs[!idx.dup],sites[!idx.dup,])
  cat("# of sites with unique IDs", nrow(out.0), "\n")

  sites.dup <- sites[idx.dup,]
  sites.dup.siteIDs <- siteIDs[idx.dup]

  sf.fun.1 <- function(idx.sel) {
    sites.sel <- sites.dup[idx.sel,]
    out.sub <- by(sites.sel,
                  list(chr.ID=paste(sites.sel$chr,
                         sites.sel$pattern,
                         sites.sel$ID,                     
                         sep="-")),
                  function(x) {
##                    nr <- dim(x)[1]
##                    if(nr<2) x
##                    else x[which.min(x[,p.value.col]),]
                    x[which.min(x[,p.value.col]),]
                  })
    ldply(out.sub, function(x) x)  
  }
 
  sites.dup.cut <- cut(seq(length(dupIDs)), n.para)
  sites.dup.cut.levels <- levels(sites.dup.cut)

  sf.fun.2 <- function(cutlevel) {
    this.idx <- which(sites.dup.cut == cutlevel)
    next.idx <- which(sites.dup.siteIDs%in%dupIDs[this.idx])
    sf.fun.1(next.idx)
  }

  out.2.list <- sfLapply(sites.dup.cut.levels, sf.fun.2)
  out.2 <- ldply(out.2.list, function(x) x)
  out <- rbind(out.0,out.2)
  cat("# of sites returned:",nrow(out),"\n")
  out
}

###. Non-sharp issues
##____________________________________________________________________
choose.site.pattern <- function(e.res, start, end, e.cutoff=-1, site.cutoff=.66) {
  ##browser()
  site.patt.df <- subset(e.res,
                         subset=as.integer(row.names(e.res))>=as.integer(start) & as.integer(row.names(e.res)) <= as.integer(end)) > e.cutoff
  patt.counts <- colSums(site.patt.df)
  n.patt <- nrow(site.patt.df)
  names.patt <- dimnames(site.patt.df)[[2]]
  patts.exist <- names.patt[patt.counts>0]
  patt.best <- which.max(patt.counts)
  patt.best.name <- names.patt[patt.best]
  patt.cutoff <- floor(n.patt*site.cutoff)
  patt.best.valid <- patt.counts[patt.best] > patt.cutoff
  out.list <- list(counts=patt.counts,
                   N=n.patt,
                   patterns=paste(patts.exist, collapse=","),
                   best.pattern=patt.best.name,
                   best.cutoff=patt.cutoff,
                   best.valid=patt.best.valid)
  out.list
}


choose.sitePattern.byChr <- function(chrom, sites, e.cutoff=-1, site.cutoff=0.66) {
  ##browser()
  get.ws(ws="pattRecog",
         what="e.res",
         chr=chrom)
  .sites <- subset(sites, subset=chr==chrom)
  ## .out <- apply(.sites, 1, function(x) {cat(x['.id'],"\t");
  ##                                       cat(x['start'], x['end'], '\n');
  ##                                       choose.site.pattern(e.res,x['start'], x['end'],
  ##                                                           e.cutoff,
  ##                                                           site.cutoff)})
  .out <- apply(.sites, 1, function(x) {choose.site.pattern(e.res, x['start'], x['end'],
                                                            e.cutoff,
                                                            site.cutoff)})

  data.frame(ldply(.out, function(x) x$counts),
             N=ldply(.out, function(x) x$N)[,2],
             best.cutoff=ldply(.out, function(x) x$best.cutoff)[,2],
             best.valid=ldply(.out, function(x) x$best.valid)[,2],
             best.pattern=ldply(.out, function(x) x$best.pattern)[,2])
}

choose.sitePattern <- function(sites, e.cutoff=-1, site.cutoff=0.7) {
  ##browser()
  chrs <- as.character(unique(sites$chr))
    ##print(chrs[i])
  .out <- sfSapply(chrs,
                   choose.sitePattern.byChr,
                   sites=sites,
                   e.cutoff=e.cutoff,
                   site.cutoff=site.cutoff,
                   simplify=F)
  ##ldply(.out, function(x) x)
  .out
}

rangedSites <- function(sites, winsize) {
  
  sites$Win <- as.integer(sites$Win)  
  sites.win <- RangesList(lapply(split(sites, sites$pattern),
                                 function(x) IRanges(start=x$Win, width=winsize)))

  s.win <- sites$Win
  s.pattern <- sites$pattern
  s.id <- sites$ID
  s.min <- tapply(s.win, list(s.id,s.pattern), min)
  s.max <- tapply(s.win, list(s.id,s.pattern), max)
  s.ends <- sapply(seq(ncol(s.min)),
                   function(icol,x,y) {idx <- is.na(x[,icol]);
                                       cbind(x[!idx,icol], y[!idx,icol]) },
                   x=s.min,
                   y=s.max)

  sites.startend <- lapply(s.ends, function(x) cbind(x[,1],x[,2]+winsize-1))
  names(sites.startend) <- dimnames(s.min)[[2]]
  
  sites.site <- RangesList(lapply(sites.startend,
                                  function(x) IRanges(start=x[,1],end=x[,2])))
  list(sites.win=sites.win,
       sites.site=sites.site)
}

site.gaps <- function(rs) {
  
  gaps.win <- with(rs,
                   gaps(sites.win,
                        start=min(start(sites.win)),
                        end=max(end(sites.win))))
  gaps.site <- with(rs, gaps(sites.site,
                             start=min(start(sites.site)),
                             end=max(end(sites.site))))
  
  gaps.within <- with(rs, intersect(sites.site, gaps.win))
  
  list(between=gaps.site,
       within=gaps.within)
}


## shear.site <- function(sites, winsize) {
##   rsites <- rangedSites(sites, winsize)
##   gaps <- site.gaps(rsites)
## }

combine.2rl <- function(x,y) {
  ## x and y are RangesList
  x.names <- names(x)
  y.names <- names(y)
  union.names <- union(x.names, y.names)
  com.names <- intersect(x.names,y.names)
  x.only.names <- setdiff(x.names, com.names)
  y.only.names <- setdiff(y.names, com.names)
  .outlist <- vector("list",length(union.names))
  names(.outlist) <- union.names
  if(length(com.names) >0) {
    for (nm in com.names) {
      .outlist[[nm]] <- c(x[[nm]],y[[nm]])
    }
  }
  com.rdl <- RangesList(.outlist)
  c(com.rdl, x[x.only.names], y[y.only.names])
}

disjoin.sites <- function(site.rl,
                          gaps.within.rl,
                          type=c("any", "start", "end", "two.ends","within","equal")) {
  type <- match.arg(type)
  if(type =="two.ends") {
    mtch1 <- subsetByOverlaps(gaps.within.rl, site.rl, type="start")
    mtch2 <- subsetByOverlaps(gaps.within.rl, site.rl, type="end")
    gaps.mtched <- combine.2rl(mtch1, mtch2)
  } else {
    gaps.mtched <- subsetByOverlaps(gaps.within.rl, site.rl, type=type)
  }
  .disjoin <- disjoin(combine.2rl(site.rl, gaps.mtched))
  out <- setdiff(.disjoin, gaps.mtched)
  out
}

patt.summary <- function(a.values, site.cutoff=0.66) {
  ##browser()
  ##a.values <- as.data.frame(values(a.rd))[,-1]
  patt.counts <- colSums(a.values)
  patt.counts.ord.12 <- sort(patt.counts, decreasing=T)[1:2]
  patt.diff.score <- abs(diff(patt.counts.ord.12))*2/sum(patt.counts.ord.12)
  n.patt <- nrow(a.values)
  names.patt <- dimnames(a.values)[[2]]
  patts.exist <- names.patt[patt.counts>0]
  patt.best <- which.max(patt.counts)
  patt.best.name <- names.patt[patt.best]
  patt.cutoff <- floor(n.patt*site.cutoff)
  patt.best.valid <- patt.counts[patt.best] > patt.cutoff
  out.list <- list(counts=patt.counts,
                   N=n.patt,
                   patterns=paste(patts.exist, collapse=","),
                   best.pattern=patt.best.name,
                   best.cutoff=patt.cutoff,
                   best.valid=patt.best.valid,
                   diff.score=patt.diff.score)
  out.list
}

pattSummary.2dataframe <- function(.out) {
  if(length(.out) ==1) {
    out <- data.frame(rbind(.out[[1]]$counts),
                      N=.out[[1]]$N,
                      best.cutoff=.out[[1]]$best.cutoff,
                      best.valid=.out[[1]]$best.valid,
                      best.pattern=.out[[1]]$best.pattern,
                      diff.score=.out[[1]]$diff.score)
                      
  } else {
    out <- data.frame(ldply(.out, function(x) x$counts)[,-1],
                    N=ldply(.out, function(x) x$N)[,-1],
                    best.cutoff=ldply(.out, function(x) x$best.cutoff)[,-1],
                    best.valid=ldply(.out, function(x) x$best.valid)[,-1],
                    best.pattern=ldply(.out, function(x) x$best.pattern)[,-1],
                    diff.score=ldply(.out, function(x) x$diff.score)[,-1])
  }
  names(out) <- c(names(out)[-(ncol(out)-(0:4))],
                  c("N","best.cutoff","best.valid","best.pattern","diff.score"))
  out
}

ranged.e.scores <- function(e.res, winsize) {
  RangedData(ranges=IRanges(start=as.integer(rownames(e.res)), width=winsize),
             e.res, space=NULL)
}

ranged.e.TF <- function(e.res, winsize, cutoff=-1) {
 RangedData(ranges=IRanges(start=as.integer(rownames(e.res)), width=winsize),
            e.res > cutoff)
}

clean.sitePatt <- function(patt, sites, e.TF, cutoff=0.5) {
  ### sites is IRanges
  ### e.TF is RangedData from ranged.e.TF
  ##browser()
  sites <- sites[[patt]]
  e.sub <- e.TF[ranges(e.TF)[[1]] %in% sites,]
  .fac <- rep(NA, length(e.sub[[1]]))
  e.ranges <- ranges(e.sub)[[1]]
  .fac <- match(e.ranges, sites)
  ## useful but not very efficient in this case.
  ##e.TF.rd <- RangedData(ranges=e.ranges, values(e.sub), space=.fac)
  ##.RDPara <- RDApplyParams(e.TF.rd, patt.summary)
  ##.summ <- rdapply(.RDPara)
  .summary <-  by(as.data.frame(values(e.sub))[,-1],
                  list(.fac),
                  function(x, cutoff) patt.summary(x, cutoff),
                  cutoff=cutoff,
                  simplify=F)
  .summary.df <- pattSummary.2dataframe(.summary)
  idx.included <- with(.summary.df,
                       best.pattern==patt & best.valid)
  sites.rd <- RangedData(sites, score=.summary.df$diff.score*log(.summary.df$N))
  ## weighted by width to avoid small site noise
  sites.rd[idx.included,]
}

compete.patt <- function(sitelist) {
  ##browser()
  olp <- findOverlaps(sitelist[[1]], sitelist[[2]])
  if(length(olp) >1) stop("Requre RangesMatchingList of length 1 between two sites")
  mm <- as.matrix(olp)
  if(nrow(mm)>0) {
    ## valid: scores in [[1]] > [[2]]
    s1 <- score(sitelist[[1]])[mm[,1]]
    s2 <- score(sitelist[[2]])[mm[,2]]
    r1 <- s1>s2
    if(sum(r1)==length(r1)) {
      sitelist[[1]] <- sitelist[[1]]
      sitelist[[2]] <- sitelist[[2]][-(mm[,2]),]
    } else if(sum(!r1)==length(r1)) {
      sitelist[[1]] <- sitelist[[1]][-(mm[,1]),]
      sitelist[[2]] <- sitelist[[2]]
    } else {
      sitelist[[1]] <- sitelist[[1]][-(mm[,1][!r1]),]
      sitelist[[2]] <- sitelist[[2]][-(mm[,2][r1]),]
    }
  }
  sitelist
}
common.expandIdx <- function(n,m,self.included=F)
  ## Function to calculate the unique combination of 1:n & 1:m.
  ## so it will list all (1,1), (1,2), (1,3), ...
  ##  but each one is a unique combination ( (1,2) and (2,1) are duplicate combination)
{
    comb <- expand.grid(1:n,1:m)
    if(self.included) comb.unique1 <- comb[comb[,1]>=comb[,2],]
    else comb.unique1 <- comb[comb[,1]>comb[,2],]
    comb.unique2 <- comb[comb[,2]>n,]
    mx <- as.matrix(rbind(comb.unique1,comb.unique2))
    mx[,rev(seq(ncol(mx)))]
}

clean.nonuniqPatt <- function(sitelist) {
  ##browser()
  l.n <- length(sitelist)
  idx <- common.expandIdx(l.n, l.n, self.included=F)
  l.nms <- names(sitelist)
  for(i in seq(nrow(idx))) {
    sitelist[idx[i,]] <- compete.patt(sitelist[idx[i,]])
    cat("comparing",l.nms[idx[i,]],"\n")
  }
  sitelist
}
  

select.siteWins <- function(patt, wins, sites, chr) {
  .wins <- wins[wins%in%sites[[patt]]]
  .fac <- match(.wins, sites[[patt]])
  .dt <- data.frame(chr, patt, cbind(.fac, start(.wins)))
  names(.dt) <- c("chr", "pattern", "ID", "Win")
  .dt
}


format.sites4diffTest <- function(wins, sites, chr) {
  .list <- lapply(names(sites[sapply(sites, length)>0]),
                  select.siteWins,
                  wins=wins,
                  sites=sites[sapply(sites, length)>0],
                  chr=chr)
  names(.list) <- names(sites[sapply(sites, length)>0])
  .out <- ldply(.list, function(x) x)[,-1]
  .out
}
## Filename:
## By: Yaomin Xu
## Date:
## Notes: BF-based methods on sites identified with DPT approach
## =============================================================================

mk.sites.report <- function(sites,
                            index.col,
                            win.col="Win",
                            win.score="Win.score",
                            score.fun=c("max","mean"),
                            winsize=50,
                            ...)
{
  rangeFun <- function(x) {
    r.x <- range(as.integer(x[,win.col]))
    s.x <- sapply(score.fun, function(y) do.call(y, list(x[, win.score],...)))
    c(r.x[1], r.x[2]+winsize, diff(r.x)+winsize, s.x)
  }
  
  ret <- ddply(sites, index.col, rangeFun)
  s.id <- apply(ret[,index.col], 1, paste, collapse="-", sep="")
  ret <- cbind(ret, s.id)
  
  scorenames <- if(length(score.fun) >1) {
    c("score", paste("score", seq(length(score.fun)-1), sep="."))
    } else "score"
  names(ret) <- c(index.col, c("start","end","size", scorenames, "siteID"))
  ret
}

mean.qint <- function(x, probs=c(1,0.75)) {
  q12 <- quantile(x, probs)
  weights <- x<=q12[1] & x>q12[2]
  weighted.mean(x, weights)
}

## e.res.null <- function(e.res, reads.poi,sites.js) {
##   browser()
##   out <- vector("list", length=length(e.res))
##   patts <- names(e.res)
##   names(out) <- patts
##   for(patt in patts) {
##     i.sites <- subset(sites.js, subset=pattern==patt)
##     if(nrow(i.sites)>0) {
##       i.TF <- !(reads.poi%in%as.numeric(i.sites$Win))
##       i.e.res.null <- e.res[[patt]][i.TF]
##       out[[patt]] <- i.e.res.null
##     } else {
##       out[[patt]] <- e.res[[patt]]
##     }
##   }
##   out
## }
    


e.res.BF.ecdf <- function(e.res, max.winsize=200, sample.size=10000, summary.fun=mean,
                          max.winsize.asym=50, sample.size.asym=1000) {
  ##browser()
  ## sample.fun <- function(ae.res,rep.n) {
  ##   ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
  ## }
  sample.fun <- function(ae.res,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ae.res, min(length(ae.res),sample.size))), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size.asym, mean=mean(ae.res), sd=sqrt(var(ae.res)/rep.n)))
    }
  }
  sapply(e.res, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
}

score.ecdf <- function(score, max.winsize=200, sample.size=10000, summary.fun=mean,
                          max.winsize.asym=50, sample.size.asym=1000) {
  ##browser()
  ## sample.fun <- function(ae.res,rep.n) {
  ##   ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
  ## }
  sample.fun <- function(ascore,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ascore, min(length(ascore),sample.size))), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size.asym, mean=mean(ascore), sd=sqrt(var(ascore)/rep.n)))
    }
  }
  sapply(score, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
}

attach.pmeans.score <- function(sites.js.bfcut.1,
                                index.col,
                                sites.js.ex,
                                pmeans.norm,
                                sample.select) {
  res1 <- pmeans.norm[sites.js.ex$Win,sample.select, drop=F]
  res2 <- data.frame(sites.js.ex, res1)
  names(res2) <- c(names(sites.js.ex), "pmeans")
  res4 <- ddply(res2, index.col, function(x) mean(x$pmeans))$V1
  res5 <- cbind(sites.js.bfcut.1, res4)
  names(res5) <- c(names(sites.js.bfcut.1), "pmeans")
  res5
}

attach.pmeans.cutoff <- function(sites.df,
                                 score.ecdf,
                                 by.chr=F,
                                 pattern.col="pattern",
                                 size.col="size",
                                 score.col="pmeans",
                                 chr.col="chr",
                                 winsize=50)
{
  ##browser()
  patts <- unique(sites.df[,pattern.col])
  chrs <- unique(sites.df[, chr.col])
  out.site <- NULL
  if(by.chr) {

  } else {
    for(ipatt in patts) {
      sites.df.patt <- subset(sites.df, subset=pattern==ipatt)
      bf.ecdf.patt <- score.ecdf[[1]]
      bfP <- apply(sites.df.patt, 1, function(x) {
        ##print(x);
        this.ecdf <- bf.ecdf.patt[[as.integer(as.numeric(x[size.col])/winsize)]];
        1-this.ecdf(as.numeric(x[score.col]))})
      iout <- cbind(sites.df.patt,bfP)
      names(iout) <- c(names(sites.df.patt), "score.P")
      out.site <- rbind(out.site, iout)
    }
  }
  out.site
}
 
  

bf.filter <- function(e.res,
                      sites.js.ex,
                      winsize,
                      bf.cutoff,
                      reads.poi,
                      pmeans.norm,
                      events.wins,
                      sample.select,
                      index.col=c("chr","pattern","ID"),
                      summary.fun=mean, sample.size=10000,
                      which.score=c("score.1","score","score.2","score.P","qvalue")) {
  ##browser()
  ## score.1 seems to be the most reasonable option
  if(missing(sample.select)) sample.select <- rep(T, ncol(pmeans.norm))
  which.score <- match.arg(which.score)
  sites.js.bfcut.1 <- mk.sites.report(sites.js.ex,
                                    index.col,
                                    "Win","Win.score",c("max","mean","min"),
                                    winsize)
  
  ##e.res.null <- e.res.null(e.res, reads.poi,sites.js.ex)

  if(all(dim(events.wins)==1)) {
    null.idx <- setdiff(row.names(pmeans.norm), sites.js.ex$Win)
    null.pmeans <- pmeans.norm[null.idx, sample.select,drop=F]
    null.pmeans.ecdf <- score.ecdf(null.pmeans, max.winsize=max(sites.js.bfcut.1$size)/winsize+1)
    sites.js.bfcut.1a <- attach.pmeans.score(sites.js.bfcut.1,
                                             index.col,
                                             sites.js.ex,
                                             pmeans.norm,
                                             sample.select)
    sites.js.bfcut.2 <- attach.pmeans.cutoff(sites.js.bfcut.1a,
                                             null.pmeans.ecdf)
    
  } else {
    sites.js.bfcut.2 <- sites.js.bfcut.1
  }
  
  ##bf.ecdf <- e.res.BF.ecdf(e.res, max.winsize=max(sites.js.bfcut.1$size)/winsize+1)

  ##sites.js.bfcut.2 <- attachP.BF.cutoff(sites.js.bfcut.1, bf.ecdf, winsize=winsize)
  ## if(which.score %in% c("score.1","score","score.2")) {
  ##   sites.BF <- subset(sites.js.bfcut.2, subset=sites.js.bfcut.2[,which.score] > bf.cutoff)
  ## } else if(which.score == "score.P") {
  ##   sites.BF <- subset(sites.js.bfcut.2, subset=sites.js.bfcut.2[,which.score] < bf.cutoff)
  ## } else stop("Invalid bf score!")
  sites.BF <- bf.filter.selSites(sites.js.bfcut.2, bf.cutoff, which.score)

  sites.ID <- apply(sites.js.ex[,c("chr","pattern","ID")],1, paste, collapse="-")

  list(sites.js.ex = sites.js.ex[sites.ID%in%sites.BF$siteID,],
       sites.bf = sites.js.bfcut.2)
}

bf.filter.selSites <- function(sites.bf, cutoff,
                               which.score=c("score.1","score","score.2","score.P","qvalue")) {
  ##browser()
  which.score <- match.arg(which.score) 
  if(which.score %in% c("score.1","score","score.2")) {
    sites <- subset(sites.bf, subset=sites.bf[,which.score] > cutoff)
  } else if(which.score == "score.P") {
    sites <- subset(sites.bf, subset=sites.bf[,which.score] < cutoff)
  } else if(which.score == "qvalue") {
    pvalues <- sites.bf[,"score.P"]
    sites.qvalue <- rep(NA, length=length(pvalues))
    p.finite <- is.finite(pvalues)
    sites.qvalue[p.finite] <- compute.qvalues(pvalues[p.finite],lambda=0)
    sites <- data.frame(sites.bf, qvalue=sites.qvalue)
    sites <- subset(sites, subset= qvalue< cutoff)
  } else stop("Invalid bf score!")
  sites
}  

BF.ecdf<- function(max.winsize=200, sample.size=10000, summary.fun=mean, by.chr=F, max.winsize.asym=50) {
  ##browser()
  sample.fun <- function(ae.res,rep.n) {
    if(rep.n < max.winsize.asym) {
      ecdf(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun))
    } else {
      ecdf(rnorm(sample.size, mean=mean(ae.res), sd=sqrt(var(ae.res)/rep.n)))
    }
  }
  chrs <- get.ws.chrs()
  winsize <- get.winsize()
  out.list <- list()
  if(by.chr) {
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      e.res.nrow <- nrow(e.res)
      if(sample.size>e.res.nrow) sample.size <- e.res.nrow
      cat(chr,":",e.res.nrow)
      sfExport("e.res","max.winsize","sample.size","summary.fun","sample.fun")
      out <- sapply(e.res, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
      names(out) <- names(e.res)
      out.list <- c(out.list, list(out))
      cat("\tDONE\n")
    }
    names(out.list) <- chrs
  } else {
    ##browser()
    all.eres <- NULL
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      for(i in seq(e.res)) {
        all.eres[[i]] <- c(all.eres[[i]], e.res[[i]])  
      }
      cat(length(all.eres[[1]]),'\n')
      if(length(all.eres[[1]])>50000) break
      }
    sfExport("all.eres","max.winsize","sample.size","summary.fun","sample.fun")
    out.list <- sapply(all.eres, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)), simplify=F)
    names(out.list) <- names(e.res)
  }
    
  out.list
}

BF.cutoff <- function(max.winsize=200, sample.size=10000, summary.fun=mean, cutoff=0.001, by.chr=F) {
  ##browser()
  sample.fun <- function(ae.res,rep.n) {
    quantile(apply(replicate(rep.n, sample(ae.res, sample.size)), 1, summary.fun), probs=1-cutoff)
  }
  chrs <- get.ws.chrs()
  winsize <- get.winsize()
  out.list <- list()
  if(by.chr) {
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      e.res.nrow <- nrow(e.res)
      if(sample.size>e.res.nrow) sample.size <- e.res.nrow
      cat(chr,":",e.res.nrow)
      sfExport("e.res","max.winsize","sample.size","summary.fun","cutoff","sample.fun")
      ##sfLibrary(snowfall)
      out <- lapply(e.res, function(x) {##sfExport("x","max.winsize","sample.size","summary.fun","cutoff","sample.fun");
        sfSapply(seq(max.winsize), function(y) sample.fun(x,y))})
      names(out) <- names(e.res)
      out.list <- c(out.list, list(out))
      cat("\tDONE\n")
    }
    names(out.list) <- chrs
  } else {
    ##browser()
    all.eres <- NULL
    for (ichr in chrs) {
      assign("chr", ichr, envir = .GlobalEnv)
      get.ws("pattRecog", "e.res")
      for(i in seq(e.res)) {
        all.eres[[i]] <- c(all.eres[[i]], e.res[[i]])  
      }
      cat(length(all.eres[[1]]),'\n')
      if(length(all.eres[[1]])>50000) break
      }
    sfExport("all.eres","max.winsize","sample.size","summary.fun","cutoff","sample.fun")
    out.list <- lapply(all.eres, function(x) sfSapply(seq(max.winsize), function(y) sample.fun(x,y)))
    ##browser()
    names(out.list) <- names(e.res)
  }
  out.1 <- sapply(out.list,
                  function(x) ldply(x, function(y) data.frame(y)),
                  simplify=F)
  out.df <- ldply(out.1,
                  function(z) {names(z) <- c('pattern','x');
                               z}
                  )
  if(by.chr) names(out.df) <- c("chr","pattern","cutoff")
  else {out.df <- out.df[,c(1,3)]
        names(out.df) <- c("pattern", "cutoff")
      }
  out.df
}



select.sites.BF.cutoff <- function(sites.df,
                                   bf.cutoff,
                                   pattern.col="pattern",
                                   size.col="size",
                                   score.col="score",
                                   chr.col="chr",
                                   winsize=50
                                   )
{
  ## Assumption: sites.df
  ## > head(sites.df)
  ##     .id   chr pattern     ID    start      end size    score  score.1
  ## 1 chr10 chr10      P1  1.100 31360650 31361400  750 15.19609 12.31469
  ## 2 chr10 chr10      P1 1.1000 14603400 14603450   50 14.20799 14.20799
  ## 3 chr10 chr10      P1 1.1001 14625350 14625400   50 14.11777 14.11777
  ## 4 chr10 chr10      P1 1.1002 14646250 14646300   50 14.30350 14.30350
  ## 5 chr10 chr10      P1 1.1003 14681400 14681450   50 14.48476 14.48476
  ## 6 chr10 chr10      P1 1.1004 14730350 14730400   50 13.88597 13.88597
  ##            siteID
  ## 1  chr10-P1-1.100
  ## 2 chr10-P1-1.1000
  ## 3 chr10-P1-1.1001
  ## 4 chr10-P1-1.1002
  ## 5 chr10-P1-1.1003
  ## 6 chr10-P1-1.1004
  ##
  ##browser()
  patts <- unique(sites.df[,pattern.col])
  chrs <- unique(sites.df[, chr.col])
  out.sites <- NULL
  if("chr" %in% names(bf.cutoff)) {
    for (ichr in chrs) {
      for (ipatt in patts) {
        isites <- subset(sites.df, subset=sites.df[,pattern.col]==ipatt & sites.df[, chr.col]==ichr)
        ibfcutoff <- subset(bf.cutoff, subset=chr==ichr & pattern == ipatt)
        out.isites <- subset(isites, subset=isites[,score.col]>ibfcutoff[isites[,size.col]/50, 'cutoff'])
        out.sites <- rbind(out.sites, out.isites)
      }
    }
  } else {
    for (ipatt in patts) {
      isites <- subset(sites.df, subset=sites.df[,pattern.col]==ipatt)
      ibfcutoff <- subset(bf.cutoff, subset=pattern == ipatt)
      out.isites <- subset(isites, subset=isites[,score.col]>ibfcutoff[isites[,size.col]/50,"cutoff"])
      out.sites <- rbind(out.sites, out.isites)
    }
  }
  out.sites
}


attachP.BF.cutoff <- function(sites.df,
                              bf.ecdf,
                              by.chr=F,
                              pattern.col="pattern",
                              size.col="size",
                              score.col="score",
                              chr.col="chr",
                              winsize=50)
{
  ##browser()
  patts <- unique(sites.df[,pattern.col])
  chrs <- unique(sites.df[, chr.col])
  out.site <- NULL
  if(by.chr) {

  } else {
    for(ipatt in patts) {
      sites.df.patt <- subset(sites.df, subset=pattern==ipatt)
      bf.ecdf.patt <- bf.ecdf[[ipatt]]
      bfP <- apply(sites.df.patt, 1, function(x) {
        ##print(x);
        this.ecdf <- bf.ecdf.patt[[as.integer(as.numeric(x[size.col])/winsize)]];
        1-this.ecdf(as.numeric(x[score.col]))})
      iout <- cbind(sites.df.patt,bfP)
      names(iout) <- c(names(sites.df.patt), "score.P")
      out.site <- rbind(out.site, iout)
    }
  }
  out.site
}
## Filename:
## By: Yaomin Xu
## Date: 2011-02-07
## Notes: Sequence/motif-based methods to sites identified by the DPT approach
## =============================================================================

get.sequence <- function(sites, ref, winsize,ext=500) {
  ### sites -  is the general sites format produced at different stages of
  ###  DPT analysis. It requires to have columns: chr, pattern, ID,
  ###  start, end.
  ### ref - reference genome.
  outseq.names <- apply(sites, 1, function(x) paste(x[c("chr", "pattern", "ID")],
                                               collapse="-", sep=""))
  outseq <- subseq(ref, start=sites[1,'start']-ext, end=sites[1,'end']+winsize-1+ext)
  for(i in seq(nrow(sites))[-1]) {
    outseq <- c(outseq, subseq(ref, start=sites[i,'start']-ext, end=sites[i,'end']+winsize-1+ext))
  }
  names(outseq) <- outseq.names
  outseq
}

get.sequence.chr <- function(sites, ref, wisize,ext=500) {
  ### Assumption: sites and ref are refering to the same chromosome.
  n.site <- nrow(sites); stopifnot(n.site>0)
  outseq <- subseq(ref,
                   start=sites[1,'start']-ext,
                   end=sites[1,'end']+winsize-1+ext)
  for(i in seq(nrow(sites))[-1]) {
    outseq <- c(outseq, subseq(ref,
                               start=sites[i,'start']-ext,
                               end=sites[i,'end']+winsize-1+ext))
  }
  ## outseq <- apply(sites, 1, function(x) subseq(ref,
  ##                                              start=as.integer(x[4]),
  ##                                              end=as.integer(x[5])-1))
  names(outseq) <- sites$siteID
  outseq
}

read.refgenome <- function(chrs, path="./FASTA/", ext=".fa") {
  cat("\n")
  reflist <- sapply(chrs, function(x) {cat(x);
                                       tst <- read.DNAStringSet(paste(path,x,ext,sep=""));
                                       cat("\tDONE\n");
                                       tst})
  names(reflist) <- chrs
  reflist
}

write.sites.seqs <- function(siteseqs, filepath) {
  stopifnot(is(siteseqs, 'list'))
  chrs <- names(siteseqs)
  if(length(siteseqs) >0) {
    cat("\n", chrs[1])
    write.XStringSet(siteseqs[[1]], filepath)
    cat("\tDONE\n")
  }
  if(length(siteseqs)>1) {
    for (i in seq(along = siteseqs)[-1]) {
      cat(chrs[i])
      write.XStringSet(siteseqs[[i]], filepath, append=T)
      cat("\tDONE\n")
    }
  }
}

attach.TsiteOLP <- function(sites.chr, Tsites) {
  
  stopifnot(length(unique(sites.chr$chr)) ==1)
  chrom <- as.character(unique(sites.chr$chr))
  cat(chrom, "\n")
  ##browser()
  overlap <- vector("logical", nrow(sites.chr))
  whichtsite <- overlapProp.tsite <- vector("numeric", nrow(sites.chr))
  overlapProp.site<- overlapSize <- vector("numeric", nrow(sites.chr))
  sites <- subset(sites.chr, select=c("start","end"))
  tsites <- subset(Tsites, subset=chr==chrom, select=c("start","end"))
  if(nrow(tsites) > 1) {
    olp <- any.overlap.sites(tsites, sites, e=0)
    ##olp$where <- as.integer(as.character(olp$where))
    olp.where <- as.numeric(unlist(strsplit(paste(as.character(olp[olp$overlap,"where"]),
                                                  collapse=",",
                                                  sep=""),
                                            split=","))
                            )

    olp.noNA <- subset(olp, subset=overlap)
    new.olp.noNA <- ldply(apply(cbind(tsite=row.names(olp.noNA),olp.noNA),
                                1,
                                function(x) data.frame(cbind(x[1],
                                                             as.numeric(unlist(strsplit(x[4], split=",")))
                                                             )
                                                       )
                                ),
                          function(y) y
                          )
    new.olp.noNA$X1 <- as.numeric(as.character(new.olp.noNA$X1))
    new.olp.noNA$X2 <- as.numeric(as.character(new.olp.noNA$X2))
    
    tsite.olp <- tsites[new.olp.noNA$X1,]
    ##sites.olp <- sites[olp[olp$overlap, "where"],]
    sites.olp <- sites[new.olp.noNA$X2,]
    ##browser()
    
    ##overlap[olp.noNA$where] <- T
    overlap[new.olp.noNA$X2] <- T
    ##overlap[olp.where] <- T
    
    ##whichtsite[olp.noNA$where] <- row.names(tsites)[olp$overlap]
    whichtsite[new.olp.noNA$X2] <- row.names(tsites)[new.olp.noNA$X1]
    olp.prop <- overlap.proportion(tsite.olp, sites.olp)
    ##overlapProp.tsite[olp.noNA$where] <- olp.prop$p1
    overlapProp.tsite[new.olp.noNA$X2] <- olp.prop$p1
    ##overlapProp.site[olp.noNA$where] <- olp.prop$p2
    overlapProp.site[new.olp.noNA$X2] <- olp.prop$p2
    overlapSize[new.olp.noNA$X2] <- overlap.size(tsite.olp, sites.olp)
  }
  out <- cbind(sites.chr,
               overlap=overlap,
               whichTsite=whichtsite,
               overlapProp.Tsite=overlapProp.tsite,
               overlapProp.site=overlapProp.site,
               overlapSize=overlapSize)
  out
}

win.movingSum <- function(ares, win.size=50, by=50) {
  pulse.n <- length(ares)
  win.n <- pulse.n - win.size + 1
  win.seq <- seq(1,win.n,by=by)
  res <- sapply(win.seq,
                function(x, dat, win.size) sum(dat[x - 1 + seq(win.size)]),
                dat = ares,
                win.size=win.size)
  res.sum <- cumsum(res)
  data.frame(winseq=win.seq+win.size-1, n=res.sum, pct=res.sum/(win.seq+win.size-1))
}

merge.fimo.site <- function(f, s) {
  
  s.1 <- subset(s, subset=hf.fimo)
  s.2 <- subset(s, subset=!hf.fimo)
  row.names(s.1) <- as.character(s.1$siteID)
  sf.1 <- s.1[as.character(f$V2),]
  full.1 <- cbind(sf.1,f)
  full.1.names <- names(full.1)

  sf.2 <- matrix(NA, nrow=nrow(s.2), ncol=ncol(f))
  full.2 <- cbind(s.2, sf.2)
  names(full.2) <- full.1.names

  full <- rbind(full.1, full.2)
  row.names(full)<- seq(nrow(full))
  full
}

extract.info.fromQuESTPeakOutput <- function(peakNames) {
  process.name <- function(aname) {
    ##browser()
    nmfull <- unlist(strsplit(aname, split="~"))
    nmout.1 <- nmfull[c(1,2,3)]
    nmout.2 <- unlist(strsplit(nmfull[9], split="-"))
    nmout.3 <- nmfull[c(17,19)]
    nmout <- c(nmout.1, nmout.2, nmout.3)
  }
  fullInfo <- ldply(sfClusterApplyLB(peakNames, process.name), function(x) x)
  row.names(fullInfo) <- peakNames
  fullInfo <- cbind(fullInfo, peakNames)
  names(fullInfo) <- c("ID", "chr","center","start","end","qv","rank","siteID")
  class(fullInfo$center) <- class(fullInfo$start) <- class(fullInfo$end) <- "numeric"
  class(fullInfo$qv) <- class(fullInfo$rank) <- "numeric"
  fullInfo
}
extract.info.fromQuESTRegioOutput <- function(peakNames) {
  ##browser()
  process.name <- function(aname) {
    ##browser()
    nmfull <- unlist(strsplit(aname, split="~"))
    nmout.1 <- nmfull[c(1,2)]
    nmout.2 <- unlist(strsplit(nmfull[3], split="-"))
    nmout.3 <- nmfull[c(9,23,25)]
    nmout <- c(nmout.1, nmout.2, nmout.3)
  }
  fullInfo <- ldply(sfClusterApplyLB(peakNames, process.name), function(x) x)
  row.names(fullInfo) <- fullInfo[,1]
  fullInfo <- cbind(fullInfo, fullInfo[,1])
  names(fullInfo) <- c("ID", "chr","start","end","center","qv","rank","siteID")
  class(fullInfo$center) <- class(fullInfo$start) <- class(fullInfo$end) <- class(fullInfo$rank)<- "integer"
  class(fullInfo$qv) <- "numeric"
  fullInfo
}

get.questRegions <- function(fasta.file) {
  qnames <- names(read.DNAStringSet(fasta.file))
  regions <- extract.info.fromQuESTRegioOutput(qnames)
  dlply(regions, 'chr', function(x) x)
}

extend.QuESTFimoNames <- function(shortNames, fullNames) {
  ##browser()
  shortNames <- as.character(shortNames)
  process.ext <- function(sname) {
    grep(sname, fullNames)
  }
  sfExport("fullNames")
  extendedNames <- unlist(sfClusterApplyLB(shortNames, process.ext))
  fullNames[unlist(extendedNames)]
}

remove.duplicateRegion <- function(x) {
  whichmin <- which.min(x$rank)
  cbind(x[whichmin,], nrow(x))
}
## Filename: analysis.pgm.R
## By: Yaomin Xu
## Date:
## Notes:Post DPT functions
## =============================================================================
clean.dupcontr.sigsites <- function(sigsites,
                                    p.value.col)
{
  ## this is a inefficient implementation, optimize it whenever possible
  out.1 <- by(sigsites,
              list(chr.ID=paste(sigsites$chr,
                     sigsites$pattern,
                     sigsites$ID,                     
                     sep=":")),
              function(x) {
                nr <- dim(x)[1]
                if(nr<2) x
                else x[which.min(x[,p.value.col]),]
              })
  out.2 <- ldply(out.1, function(x) x)
  out.2
}

map.regs.chr <- function(trf.sig,
                         trf,
                         chrom,
                         methylType,
                         patts=c("**1","*1*","1**"),
                         cutoff=5,
                         winsize=100)
{
  
  trf.sig.chr <- subset(trf.sig, subset=chr==chrom & MethylType==methylType)
  trf.sig.ord <- trf.sig.chr[order(trf.sig.chr$start),]

  out <- trf.sig.ord[,c('ID','pattern')]
  for(i in seq(patts)) {
    patt <- patts[i]
    trf.s.part <- subset(trf,
                         subset= (pattern == patt &
                                  chr==chrom)
                         )
    trf.s <- cbind(trf.s.part, size=with(trf.s.part, end - start))
    trf.s.c <- subset(trf.s, subset = size > cutoff*winsize)
    trf.s.c.ord <- trf.s.c[order(trf.s.c$start),]

    idx <- any.overlap.sites(trf.sig.ord[,c('start','end')],
                             trf.s.c.ord[,c('start','end')],
                             e = 0,
                             self=F)
    unique.overlap <- idx$count==1
    trf.sig.idx <- seq(nrow(idx))[unique.overlap]
    ## trf.idx <- as.numeric(unlist(strsplit(paste(idx[idx$overlap,]$where,
    ##                                             collapse=","),
    ##                                       split=",")))
    trf.idx <- as.numeric(as.character(idx[unique.overlap,'where']))
    
    out.here <- matrix(rep(NA, 2*nrow(trf.sig.ord)), ncol=2)
    out.here[trf.sig.idx,1] <- trf.s.c.ord[trf.idx, 'ID']
    out.here[trf.sig.idx,2] <- as.character(trf.s.c.ord[trf.idx, 'pattern'])
    
    out <- cbind(out, out.here)
  }
  names(out) <- c('s.ID', "s.Pattern",
                  sapply(seq(patts), function(x) paste(c('r.ID','r.Pattern'),x, sep=".")))
  out
}

overlapProp.byID <- function(dat1,dat2,id1,id2) {
  ##browser()
  stopifnot(nrow(id1)==nrow(id2))
  p1 <- p2 <- rep(NA, nrow(id1))
  na1 <- is.na(id1[,1])
  na2 <- is.na(id2[,1])
  nona <- !(na1 | na2)

  ## subsetting
  ## dat1.idx <- match(paste(id1[,1],id1[,2],sep="."),with(dat1,paste(ID,pattern,sep=".")))
  ## dat2.idx <- match(paste(id2[,1],id2[,2],sep="."),with(dat2,paste(ID,pattern,sep=".")))

  ## dat1.sub <- dat1[dat1.idx,]
  ## dat2.sub <- dat2[dat2.idx,]
  dat.commSub <- overlap.commonSet(dat1, dat2, id1, id2)
  ##id1.nna <- match(id1[nona], dat1$ID)
  ##id2.nna <- match(id2[nona], dat2$ID)

  R1 <- dat.commSub$dat1[nona,c('start','end')]
  R2 <- dat.commSub$dat2[nona,c('start','end')]
  props <- overlap.propotion(R1,R2)

  p1[nona] <- props$p1
  p2[nona] <- props$p2
  data.frame(p1=p1, p2=p2)
}

overlap.commonSet <- function(dat1, dat2, id1, id2) {
  stopifnot(nrow(id1)==nrow(id2))
  ## subsetting
  dat1.idx <- match(paste(id1[,1],id1[,2],sep="."),with(dat1,paste(ID,pattern,sep=".")))
  dat2.idx <- match(paste(id2[,1],id2[,2],sep="."),with(dat2,paste(ID,pattern,sep=".")))

  dat1.sub <- dat1[dat1.idx,]
  dat2.sub <- dat2[dat2.idx,]

  list(dat1=dat1.sub, dat2=dat2.sub)
}

map.regions <- function(report.sig,
                        report,
                        chrom,
                        methylType,
                        patts=c("**1","*1*"),
                        cutoff=5,
                        winsize=100)
{  ## parallel the process on multiple processes
  trf.sig <- subset(report.sig,
                    subset= chr==chrom & MethylType==methylType)
  trf <- subset(report,
                subset= chr== chrom)

  reg.map <- map.regs.chr(report.sig,
                          report,
                          chrom,
                          methylType,
                          patts=patts,
                          cutoff=cutoff,
                          winsize=winsize)

  tst <- overlapProp.byID(trf.sig, trf, reg.map[,1:2], reg.map[,3:4])
  tst.subset <- overlap.commonSet(trf.sig, trf, reg.map[,1:2], reg.map[,3:4])

  tst1 <- overlapProp.byID(trf.sig, trf, reg.map[,1:2], reg.map[,5:6])
  tst1.subset <- overlap.commonSet(trf.sig, trf, reg.map[,1:2], reg.map[,5:6])

  tst.p2df <- data.frame(prop1=tst$p2, prop2=tst1$p2)

  tst.p2df.na <- apply(tst.p2df, 1, function(x) all(is.na(x)))

  tst.p2df.nona <- tst.p2df[!tst.p2df.na,]

  tst.p2df.choose <- rep(0, nrow(tst.p2df))

  tst.p2df.choose[!tst.p2df.na] <- apply(tst.p2df.nona, 1, which.max)


  tst.regions <- cbind(tst.subset$dat1,
                       (!tst.p2df.na),
                       tst.p2df,
                       tst.subset$dat2,
                       reg.map)

  tst.regions[tst.p2df.choose==0,] <- data.frame(
                tst.subset$dat1[tst.p2df.choose==0,],
                (!tst.p2df.na)[tst.p2df.choose==0],
                tst.p2df[tst.p2df.choose==0,],
                tst.subset$dat1[tst.p2df.choose==0,1:12],
                reg.map[tst.p2df.choose==0,]
                )
  tst.regions[tst.p2df.choose==1,] <- cbind(
                tst.subset$dat1[tst.p2df.choose==1,],
                (!tst.p2df.na)[tst.p2df.choose==1],
                tst.p2df[tst.p2df.choose==1,],
                tst.subset$dat2[tst.p2df.choose==1,],
                reg.map[tst.p2df.choose==1,])
  tst.regions[tst.p2df.choose==2,] <- cbind(
                tst.subset$dat1[tst.p2df.choose==2,],
                (!tst.p2df.na)[tst.p2df.choose==2],
                tst.p2df[tst.p2df.choose==2,],
                tst1.subset$dat2[tst.p2df.choose==2,],
                reg.map[tst.p2df.choose==2,])
  names(tst.regions) <- c(names(tst.subset$dat1),
                          "R.exist",
                          paste("R", names(tst.p2df), sep="."),
                          paste("R", names(tst.subset$dat2), sep="."),
                          names(reg.map))
  tst.regions
}
###... Straight forward approach for region identification

map.region.chr.2 <- function(tst, gap=500) {
 
  tst.o <- tst[order(as.numeric(tst$start)),]
  tst.d <- tst.o$start[-1]-tst.o$end[-nrow(tst.o)]
  tst.flw <- tst.d < gap

  tst.id <- rep(0, nrow(tst))
  for(j in seq(nrow(tst.o))) {
    if(j==1) {
      tst.id[j] <- 1
    }
    else {
      if(tst.flw[j-1]) tst.id[j] <- tst.id[j-1]
      else {
        tst.id[j] <- tst.id[j-1] + 1
      }
    }
  }

  r.id <- paste("R", tst.o$chr, tst.id, sep="-")

  transform(tst.o, R.id=r.id)
  
  
}

map.region.2 <- function(tst, gap=500) {
  chrs <- unique(tst$chr)
  out.1 <- NULL
  for(chrom in chrs) {
    tst.chr <- subset(tst, subset=chr==chrom)
    out.1 <- rbind(out.1, map.region.chr.2(tst.chr, gap))
  }
  
  r.range <- t(sapply(by(out.1[,c('start','end')], out.1$R.id, range), function(x) x))
  out.2 <- cbind(out.1,r.range[match(out.1$R.id, row.names(r.range)),])
  names(out.2) <- c(names(out.1), c('r.start','r.end'))
  out.2
}

add.siteUID <- function(tst) {
  sid <- paste(tst$chr, tst$pattern, tst$ID, sep="-")
  out <- transform(s.id=sid, tst)
  row.names(out) <- sid
  out
}

map.region <- function(test.report, gap) {
  ##browser()
### test.report has chr, start, and end columns
### gap is the max gap that will get joined, so the minimum gap after joining will be gap+1
  names(test.report) <- replace.names(names(test.report), "chr", "space")
  .rd <- as(test.report, "RangedData")
  .rl <- ranges(.rd)
  region.rl <- reduce(.rl, drop.empty.ranges=T, min.gapwidth=gap+1)
  ovlp <- as.matrix(findOverlaps(.rl, region.rl))
  r.id <- paste(.rd$space, "r",ovlp[,2], sep="_")
  .rd$r.id <- r.id
  out <- as.data.frame(.rd)
  names(out) <- replace.names(names(out), "space", "chr")
  out

}

calculate.gmeans <- function(tst.regions) {
  
  mx01 <-  ldply(strsplit(as.character(tst.regions$contrast), split=character()), function(x) x)
  n0 <- apply(mx01, 1, function(x) sum(x=='0'))
  n1 <- apply(mx01, 1, function(x) sum(x=='1'))
  n <- n0+n1
  y1 <- tst.regions$value.0+n0/n*tst.regions$value.1
  y0 <- y1 - tst.regions$value.1
  mean0 <- y0
  mean1 <- y1
  ##mean0[tst.regions$value.1>0] <- y0[tst.regions$value.1>0]
  mean0[tst.regions$value.1<=0] <- y1[tst.regions$value.1<=0]
  ##mean1[tst.regions$value.1>0] <- y1[tst.regions$value.1>0]
  mean1[tst.regions$value.1<=0] <- y0[tst.regions$value.1<=0]

  data.frame(mean.0=mean0,
             mean.1=mean1)
}

cleanup.on.gmeans <- function(regions) {
  
  cat("clean up non-sharp sites ...\n")
  regions$patt <- regions$pattern
  uniquepatts <- as.character(unique(regions$patt))
  patt111 <- uniquepatts[sapply(strsplit(uniquepatts,
                                         split=""),
                                function(x) all(x=="1"))]
  ###... remaking the 0/1 decision based on LDA
  training <- subset(regions, subset=pattern!=patt111)
  lda.m <- with(training, lda(as.matrix(c(mean.0, mean.1)),
                              grouping=c(rep(0, length(mean.0)), rep(1, length(mean.1)))))
  ###... update site patterns based on LDA model
  regions.hyper.no111 <- subset(regions, subset=MethylType=="Hyper"&pattern!=patt111)
  sharp.hyper.no111 <- with(regions.hyper.no111,
                            predict(lda.m, as.matrix(mean.0))$class==0 & predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hyper.no111$pattern[!sharp.hyper.no111] <- patt111
  cat("num of validated sharp sites in hyper not-all-1 patterns :",
      sum(sharp.hyper.no111),
      "/",
      nrow(regions.hyper.no111) ,
      "\n")
  
  regions.hyper.111 <- subset(regions, subset=MethylType=="Hyper"&pattern==patt111)
  sharp.hyper.111 <- with(regions.hyper.111,
                          predict(lda.m, as.matrix(mean.0))$class==0 & predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hyper.111$pattern[sharp.hyper.111] <- regions.hyper.111$contrast[sharp.hyper.111]
  cat("num of validated sharp sites in hyper all-1 patterns :",
      sum(sharp.hyper.111),
      "/",
      nrow(regions.hyper.111) ,
      "\n")
   
  regions.hypo.no111 <- subset(regions, subset=MethylType=="Hypo"&pattern!=patt111)
  sharp.hypo.no111 <- with(regions.hypo.no111,
                           predict(lda.m, as.matrix(mean.0))$class==0&predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hypo.no111$pattern[!sharp.hypo.no111] <- patt111
  cat("num of validated sharp sites in hypo not-all-1 patterns :",
      sum(sharp.hypo.no111),
      "/",
      nrow(regions.hypo.no111) ,
      "\n")
  

  regions.hypo.111 <- subset(regions, subset=MethylType=="Hypo"&pattern==patt111)
  sharp.hypo.111 <- with(regions.hypo.111,
                         predict(lda.m, as.matrix(mean.0))$class==0&predict(lda.m, as.matrix(mean.1))$class==1)
  regions.hypo.111.patts <- regions.hypo.111$contrast[sharp.hypo.111]
  regions.hypo.111.patterns <- gsub("b","0", gsub("0", "1", gsub("1","b", regions.hypo.111.patts)))
  regions.hypo.111$pattern[sharp.hypo.111] <- regions.hypo.111.patterns
  cat("num of validated sharp sites in hypo all-1 patterns :",
      sum(sharp.hypo.111),
      "/",
      nrow(regions.hypo.111) ,
      "\n")
  
  rbind(regions.hyper.no111, regions.hyper.111, regions.hypo.no111, regions.hypo.111)
}

get.e.res <- function(start, end, chr) {
  get.ws("pattRecog", what=c("e.res"), chr=chr)
  poi <- as.integer(rownames(e.res))
  .out <- e.res[poi>=start&poi<=end,]
  .out[order(as.integer(rownames(.out))),]
}


select.reads <- function(reads.df,
                         sample.select) {  
  poi <- reads.df$poi
  s.id <- reads.df$s.id
  cols.rmv <- match(c("poi", "s.id"), names(reads.df))
  if(missing(sample.select)) {
    reads.sub <- reads.df[,-cols.rmv]
  } else {
    reads.sub <- reads.df[,-cols.rmv][,rep(sample.select,2)]
  }
  cbind(poi=poi, reads.sub, s.id=s.id)
}

