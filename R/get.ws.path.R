#' @export
.get.ws.path <-
function(which=c("root",
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
