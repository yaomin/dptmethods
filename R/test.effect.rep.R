test.effect.rep <-
function(wins, psig, vars, group.label, rep.label,
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
