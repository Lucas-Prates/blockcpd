#include <Rcpp.h>

using namespace Rcpp;

//' Sampler for depedent models
//[[Rcpp::export]]
NumericMatrix rcpd_cpp(String family,
              int n,
              int m,
              IntegerVector  changepoints,
              List parameters){

  NumericMatrix data_mat(n, m);
  if(family == "binaryMarkov"){
    data_mat(_, 0) = Rcpp::rbinom(n, 1, 0.5);
    NumericVector p00 = parameters(0);
    NumericVector p11 = parameters(1);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < changepoints.size() - 1; j ++){
        for(int c = changepoints[j]; c < changepoints[j + 1]; c++){
          if(c > 0){
            if(data_mat(i, c - 1) == 0){
              data_mat(i, c) = R::runif(0, 1) <  p00[j] ? 0 : 1;
            }else{
              data_mat(i, c) = R::runif(0, 1) <  p11[j] ? 1 : 0;
            }
          }
        }
      }
    }
  }

  return(data_mat);
}
