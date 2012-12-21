#' @export
get.wsoption.path <-
function(which=c("root",
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
