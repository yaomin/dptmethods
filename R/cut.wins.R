cut.wins <-
function(wins.df, h, cutoff.cat=3000, cutoff.last=2000)
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
