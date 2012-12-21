match.mx <-
function(x.mx, table.mx) {
  if(is.null(dim(table.mx))) dim(table.mx) <- c(1,length(table.mx))
  reslt <- rep(0, nrow(x.mx))
  reslt <- unlist(apply(x.mx, 1, function(a.x.mx, table.mx) {
                            rst <- which(apply(table.mx, 1, 
                                        function (a.table.mx, a.x.mx) all(a.x.mx == a.table.mx),
                                        a.x.mx = a.x.mx                           
                                       )
                                );
                            rst <- ifelse(length(rst)>0, rst,NA);
                            rst
                            },
                table.mx = table.mx
               ),
               
               )
  reslt
}
