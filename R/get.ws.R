#' @export
get.ws <-
function(ws=c("mixPoisson","pattRecog","diffTest","joinSample","log","report","analysis"),
                   what,
                   chr,
                   ws.path,
                   pos=parent.frame(),
                   pos.out=parent.frame())
{
  ##output.list <- NULL
  ws <- match.arg(ws)  
  ##if(missing(pos.out)) pos.out <- parent.frame()
  if(ws == "mixPoisson") {
    if(missing(chr)) chr <- get('chr', pos=pos)
    if(missing(ws.path)) ws.path <- get.ws.path("mixPoisson", pos=pos)
    load.mixPoiData(ws.path, chr, pos=pos)
    for (eachwhat in what) {
      ##output.list <- list(output.list, get(eachwhat))
      assign(eachwhat, get(eachwhat), pos=pos.out) 
    }
  }
  else if(ws == "pattRecog") {
    if(missing(chr)) chr <- get('chr', pos=pos)
    if(missing(ws.path)) ws.path <- get.ws.path("pattRecog", pos=pos)
    load.pattRecogData(ws.path, chr, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
  }
  else if(ws == "diffTest") {
   if(missing(chr)) chr <- get('chr', pos=pos)
   if(missing(ws.path)) ws.path <- get.ws.path("diffTest", pos=pos)
   load.diffTestData(ws.path, chr, pos=pos)
   for (eachwhat in what) {
     assign(eachwhat, get(eachwhat), pos=pos.out)
   }
 }
  else if(ws == "report") {
    if(missing(ws.path)) ws.path <- get.ws.path("report", pos=pos)
    load.reportData(ws.path, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
    
  }
  else if(ws == "analysis") {
    if(missing(ws.path)) ws.path <- get.ws.path("report", pos=pos)
    load.analysisData(ws.path, pos=pos)
    for (eachwhat in what) {
      assign(eachwhat, get(eachwhat), pos=pos.out)
    }
    
  }
  ##names(output.list) <- what
  ##output.list
}
