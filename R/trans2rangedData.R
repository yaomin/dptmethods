trans2rangedData <-
function(df.in, ext.end=49) {
  ## df.in needs to have start, end, and chr as columns
  ## end in df.in is the start position of the last window
  ## ext.end set to NULL for properlly calculated ends.
  if(!is.null(ext.end)) df.in$end <- df.in$end+ext.end
  df.in <- cbind(df.in, space = df.in$chr)
  as(df.in, "RangedData")
}
