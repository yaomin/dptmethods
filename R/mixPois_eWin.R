#' @export
mixPois_eWin <-
function(y, initials=NULL, para.priors, controls,
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
