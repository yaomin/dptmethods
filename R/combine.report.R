combine.report <-
function(wins, tests, tests.select) {
  ##browser()
  if(missing(tests.select)) tests.select=grep("ID|Value|p.value", names(tests))
  part1 <- as.data.frame(t(as.data.frame(lapply(by(as.integer(wins$Win),wins$ID, range), unlist))))
  rownames.part1 <- gsub('X', '', row.names(part1))
  row.names(part1) <- rownames.part1

  tests.sel <- tests[,tests.select]
  common.rownames <- intersect(rownames.part1, row.names(tests.sel))
  part2 <- tests.sel[common.rownames,]
  
  res <- data.frame(part1[common.rownames,], part2)
  ##row.names(res) <- common.rownames
  names(res) <- c(c("start","end"), names(part2))
  res
}
