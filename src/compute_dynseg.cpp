#include <Rcpp.h>
#include "Dynseg.h"
using namespace Rcpp;


// Dynamical programming Algorithm
// Wrapper for the Dynseg class. It is an implementation of the dynamic
// programming algorithm to estimate the penalized loss estimator.
// Returns the model info
// [[Rcpp::export]]
List compute_dynseg_cpp(const List& suff_stats,
                        const String& family,
                        const int& ncol,
                        int max_blocks,
                        const Function& pen_func){

  Dynseg dynseg(family, suff_stats, pen_func, ncol, max_blocks);
  dynseg.fit_dynseg();
  List model_info(5);
  model_info[0] = dynseg.changepoints;
  model_info[1] = dynseg.parameters;
  model_info[2] = dynseg.loss; // negloglike + reg
  model_info[3] = dynseg.changepoints.size(); // number of cp
  model_info[4] = dynseg.negloglike; // negloglike

  return model_info;
}
