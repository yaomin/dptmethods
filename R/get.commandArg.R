get.commandArg <-
function(key, mode=c("numeric","character")) {
  mode <- match.arg(mode)
  if(mode == "numeric") argInput <- as.numeric(gsub(key,"", gsub('-','',grep(key, commandArgs(), v=T))))
  else argInput <- gsub(key,"", gsub('-','',grep(key, commandArgs(), v=T)))
  if(length(argInput)==0) stop("Please provide a command argument as '--cutoff-0.01'!")
  argInput
}
