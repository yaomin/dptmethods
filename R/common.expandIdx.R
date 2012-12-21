common.expandIdx <-
function(n,m,self.included=F)
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
