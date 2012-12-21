which.zeroPattern <-
function(pattNames, cap="P") {
  which(as.numeric(gsub(cap, "", pattNames))==0)
}
