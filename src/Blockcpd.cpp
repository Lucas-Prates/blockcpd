#include <Rcpp.h>
#include "Blockcpd.h"
using namespace Rcpp;

Blockcpd::Blockcpd(String family, const List& suff_stats, Function pen_func, int ncol)
  : family(family),  suff_stats(suff_stats), pen_func(pen_func), ncol(ncol),
    changepoints(0), loss(0), negloglike(0) {}

float Blockcpd::compute_negloglike(const int& left_index,
                                   const int& right_index){

  double loglike = 0;
  int block_size = 0;

  // log likelihood computation for bernoulli distribution
  if(family == "bernoulli"){
    double prob = 0;
    NumericVector col_sum = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    if(left_index - 1 != 0){
      prob = col_sum[right_index - 1] - col_sum[left_index - 2];
      block_size = col_samples[right_index - 1] - col_samples[left_index - 2];
    } else {
      prob = col_sum[right_index - 1];
      block_size = col_samples[right_index - 1];
    }
    prob /= block_size;
    if( (prob != 0)&(prob != 1) ) {
      loglike = block_size*(prob*log(prob) + (1-prob)*log(1-prob));
    }
  }

  // log likelihood computation for normal distribution
  if(family == "normal"){
    double block_mean = 0;
    double block_squared_mean = 0;
    NumericVector col_sum = suff_stats[0];
    NumericVector col_squared_sum = suff_stats[1];
    NumericVector col_samples = suff_stats[2];
    if(left_index - 1 != 0){
      block_mean = col_sum[right_index - 1] - col_sum[left_index - 2];
      block_squared_mean = col_squared_sum[right_index - 1] - col_squared_sum[left_index - 2];
      block_size = col_samples[right_index - 1] - col_samples[left_index - 2];
    } else {
      block_mean = col_sum[right_index - 1];
      block_squared_mean = col_squared_sum[right_index - 1];
      block_size = col_samples[right_index - 1];
    }
    block_mean /= block_size;
    block_squared_mean /= block_size;
    // loglike part that depends on parameters
    loglike = -block_size*log(block_squared_mean - block_mean*block_mean);
    // adding constant part to loglike
    loglike -= (block_size/2.0)*(log(2*PI) + 1);
  }

  return(-loglike);
}

float Blockcpd::compute_loss(const int& left_index,
                             const int& right_index){

  return compute_negloglike(left_index, right_index) + *REAL(pen_func(left_index, right_index));

}

void Blockcpd::fit_family_parameters(){
  int n_blocks = changepoints.size() + 1;
  IntegerVector block_begin(n_blocks), block_end(n_blocks);
  // Compute the parameters of each block estimated
  for(int i = 0; i < n_blocks; i++) {
    //Start of each block
    if(i == 0) {
      block_begin[i] = 1;
    } else {
      block_begin[i] = changepoints[i-1] + 1;
    }

    // End of each block
    if(i == n_blocks-1) {
      block_end[i] = ncol;
    } else {
      block_end[i] = changepoints[i];
    }

  }

  // Fitting parameters for bernoulli distribution
  if(family == "bernoulli"){
    NumericVector col_sums = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    int block_size;
    NumericVector probs(n_blocks);
    for(int i = 0; i < n_blocks; i++){
      probs[i] = 0;
      block_size = 0;
      if(block_begin[i] - 1 != 0){
        probs[i] = col_sums[block_end[i] - 1] - col_sums[block_begin[i] - 2];
        block_size = col_samples[block_end[i] - 1] - col_samples[block_begin[i] - 2];
      } else {
        probs[i] = col_sums[block_end[i] - 1];
        block_size = col_samples[block_end[i] - 1];
      }
      probs[i] /= block_size;
    }
    parameters = List::create(Named("probs") = probs);
  }

  // Fitting parameters for normal distribution
  if(family == "normal"){
    NumericVector block_mean(n_blocks);
    NumericVector block_var(n_blocks);
    NumericVector col_sums = suff_stats[0];
    NumericVector col_squared_sum = suff_stats[1];
    NumericVector col_samples = suff_stats[2];
    int block_size;
    for(int i = 0; i < n_blocks; i++){
      block_mean[i] = 0;
      block_var[i] = 0;
      block_size = 0;
      if(block_begin[i] - 1 != 0){
        block_mean[i] = col_sums[block_end[i] - 1] - col_sums[block_begin[i] - 2];
        block_var[i] = col_squared_sum[block_end[i] - 1] - col_squared_sum[block_begin[i] - 2];
        block_size = col_samples[block_end[i] - 1] - col_samples[block_begin[i] - 2];
      } else {
        block_mean[i] = col_sums[block_end[i] - 1];
        block_var[i] = col_squared_sum[block_end[i] - 1];
        block_size = col_samples[block_end[i] - 1];
      }
      block_mean[i] /= block_size;
      block_var[i] /= block_size;
      block_var[i] -= block_mean[i]*block_mean[i];
    }
    parameters = List::create(Named("mean") = block_mean,
                              Named("variance") = block_var);
  }

}

void Blockcpd::sort_changepoints(){
  std::sort(changepoints.begin(), changepoints.end());
}

