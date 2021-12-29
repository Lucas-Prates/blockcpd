#include <Rcpp.h>
#include "Hierseg.h"
using namespace Rcpp;


// Wrapper for the Hierarchical segmentation Algorithm cpp
//
// Obtain the solution for the optimization problem using the
// hierarchical algorithm. Creates an object of the hierseg model, fits it and
// wraps up its output
// [[Rcpp::export]]
List compute_hierseg_cpp(const List& suff_stats,
                         const String& family,
                         const int& ncol,
                         const int& min_block_size,
                         const Function& pen_func) {

  Hierseg hierseg(family, suff_stats, pen_func, ncol, min_block_size);

  // Estimation of parameters
  hierseg.fit_hierseg();

  List model_info(4);
  model_info[0] = hierseg.changepoints;
  model_info[1] = hierseg.parameters;
  model_info[2] = hierseg.loss;
  model_info[3] = hierseg.negloglike;
  return model_info;
}

