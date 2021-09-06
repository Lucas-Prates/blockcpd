#include <Rcpp.h>
using namespace Rcpp;

// Compute sufficient statistics for distributions where the R alternative
// might be slow
//[[Rcpp::export]]
List compute_suff_stats_cpp(const IntegerMatrix& data_mat,
                            const String& family){

  List suff_stats;
  if(family == "binaryMarkov") {
    int n = data_mat.nrow();
    int m = data_mat.ncol();

    // Number of samples with starting state 0 and 1,
    // and number of samples that goes from 0 to 0 and from 1 to 1
    // Initialized with zeros by default
    NumericVector N0(m - 1), N1(m - 1), p00(m - 1), p11(m - 1);

    for(int j = 0; j < m - 1; j++){
      // Cumulative sum of values
      // Case 0 is zero by default initialization
      if(j > 0){
        N0[j] = N0[j - 1];
        p00[j] = p00[j - 1];
        N1[j] = N1[j - 1];
        p11[j] = p11[j - 1];
      }
      for(int i = 0; i < n; i++){
        // When current state or the next state is NA, the iteration is skipped

        // When current state is 0
        if(data_mat(i, j) == 0){
          if(data_mat(i, j + 1) == 0){
            N0[j]++;
            p00[j]++;
          }
          if(data_mat(i, j + 1) == 1){
            N0[j]++;
          }
        }
        // When current state is 1
        if(data_mat(i, j) == 1){
          if(data_mat(i, j + 1) == 1){
            N1[j]++;
            p11[j]++;
          }
          if(data_mat(i, j + 1) == 0){
            N1[j]++;
          }
        }
      }
    }

  List bm_suff_stats(4);
  bm_suff_stats(0) = N0; bm_suff_stats(1) = p00;
  bm_suff_stats(2) = N1; bm_suff_stats(3) = p11;
  suff_stats = bm_suff_stats;
  }

  return suff_stats;
}
