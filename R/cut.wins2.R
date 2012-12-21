cut.wins2 <-
function(wins.df, h, winsize) {
  ### partially done what's in cut.wins with better performance
  ### smaller trees are not merged into the bigger branch
  ### rewrite win.search
  ### NOT READY FOR USE
  .ir <- IRanges(start=wins.df$wins, width=winsize)
  .ir.gaps <- gaps(.ir)
  .gaps.cut <- .ir.gaps[width(.ir.gaps)>h]
  split.fac <- cut(wins.df$wins,
                   breaks=sort(c(min(wins.df$wins),
                     start(.gaps.cut),
                     end(.gaps.cut),
                     max(wins.df$wins))),
                   include.lowest=T)
  split(wins.df, split.fac, drop=T)
}
