#' @export
save.ws <-
function(ws=names(dpt.options[['output']]),
                    what,
                    outpath,
                    chr=NULL,
                    pos=parent.frame()) {
  ws <- match.arg(ws)
  dpt.options <- get("dpt.options", pos=pos)
  if(missing(outpath)) outpath <- get.ws.path(ws, pos=pos)
  if(is.null(chr)) {
    outfile <- file.path(outpath, paste(dpt.options[[c('output',ws)]],"RData", sep="."))
  } else {
    outfile <- file.path(outpath, paste(dpt.options[[c('output',ws)]],chr,"RData", sep="."))
  }
  ##save.envir <- if(pos==-1) globalenv() else as.environment(pos)
  save(list=what,
       file=outfile,
       envir=as.environment(pos))
  cat(what,"from",ws,"\t","saved!\n")
}
