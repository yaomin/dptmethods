get.count.tb.sample <-
function(mappedData) {
  tb.df <- as.data.frame(table(mappedData[,"count"]))
  names(tb.df) <- c("count", "freq")
  tb.df
}
