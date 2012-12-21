Rcpp.merger <-
function(ids, smp) {
  src <-'
     string s1, s2;
     Rcpp::CharacterVector xid(id);
     Rcpp::NumericVector xpoi(poi);
     Rcpp::CharacterVector xsid(sid);
     Rcpp::NumericVector xspoi(spoi);
     Rcpp::NumericVector xscnt(scnt);
     int n_xid = xid.size();
     int n_xsid = xsid.size();
     Rcpp::NumericVector xcnt(n_xid);
     for (int i=0; i<n_xid; i++)
        for (int j=0; j<n_xsid; j++) {
          s1 = xid[i];
          s2 = xsid[j];
          if (s1 == s2 & xpoi[i]- xspoi[j] <1) xcnt[i] = xscnt[j];
          else xcnt[i] =0;
        };  
     return xcnt;
  '

  .merger <- cxxfunction(signature(id="character",
                                   poi="numeric",
                                   sid="character",
                                   spoi="numeric",
                                   scnt="numeric"),
                         includes=c("#include <string>", "using namespace std;"),
                         src,
                         plugin="Rcpp")

  .merger(ids$chr, ids$start, smp$chr, smp$start, smp$count)
}
