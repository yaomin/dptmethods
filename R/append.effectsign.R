#'  @export
append.effectsign <-
function(test.report,
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
