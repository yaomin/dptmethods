create.ws.dir <-
function(which=c("root",
                                "mixPoisson",
                                "pattRecog",
                                "diffTest",
                                "preprocess",
                                "log",
                                "report",
                                "src",
                                "tmp"),
                          pos=-1
                          )
{
  which <- match.arg(which)
  .dir <- get.wsoption.path(which, pos=pos)
  if(!file_test("-d", .dir)) {
    .create.yes <- dir.create(.dir)
    if(!.create.yes) stop(paste("Can not create",
                                .dir,
                                "for dpt analyses!",
                                sep=" "))
  } else {
    cat(.dir, "exists!","\n")
  }
}
