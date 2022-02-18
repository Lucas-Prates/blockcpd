#include <Rcpp.h>
#include "Blockcpd.h"
#include <math.h>
using namespace Rcpp;

Blockcpd::Blockcpd(String family, const List& suff_stats, Function pen_func,
                   int ncol, int min_block_size, int max_blocks)
  : family(family),  suff_stats(suff_stats), pen_func(pen_func), ncol(ncol),
    min_block_size(min_block_size), max_blocks(max_blocks),
    changepoints(0), loss(0), negloglike(0) {}

float Blockcpd::compute_negloglike(const int& left_index,
                                   const int& right_index){
  if(right_index - left_index + 1 < min_block_size){
    return INFINITY; //-loglike = -inf to force min_blocK_size restriction
  }
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
    loglike -= (block_size/2.0)*(log(2*M_PI) + 1);
  }

  if(family == "binaryMarkov"){
    if(right_index == left_index){
      return(INFINITY); // force model not to allow 1 sized blocks
    }
    else{
      NumericVector col_samples0 = suff_stats[0];
      NumericVector col_sum0 = suff_stats[1];
      NumericVector col_samples1 = suff_stats[2];
      NumericVector col_sum1 = suff_stats[3];
      double prob00 = 0, prob11 = 0;
      int n_samples0 = 0, n_samples1 = 0;
      // 1 cannot be a change point in this case!
      if(left_index - 1 > 1){
        // the -3 on the left_index - 3 is due to the convention that
        // at the end of a block, it transits to the beginning of a new block
        // using the transition matrix of the next block!!!
        prob00 = col_sum0[right_index - 2] - col_sum0[left_index - 3];
        prob11 = col_sum1[right_index - 2] - col_sum1[left_index - 3];
        n_samples0 = col_samples0[right_index - 2] - col_samples0[left_index - 3];
        n_samples1 = col_samples1[right_index - 2] - col_samples1[left_index - 3];
      } else {
        prob00 = col_sum0[right_index - 2];
        prob11 = col_sum1[right_index - 2];
        n_samples0 = col_samples0[right_index - 2];
        n_samples1 = col_samples1[right_index - 2];
      }
      if(n_samples0 > 0){prob00 /= n_samples0;}
      if(n_samples1 > 0){prob11 /= n_samples1;}

      if( (prob00 != 0)&(prob00 != 1) ) {
        loglike += n_samples0*(prob00*log(prob00) + (1-prob00)*log(1-prob00));
      }
      if( (prob11 != 0)&(prob11 != 1) ) {
        loglike += n_samples1*(prob11*log(prob11) + (1-prob11)*log(1-prob11));
      }
    }
  }
  // log likelihood computation for exponential distribution
  if(family == "exponential"){
    double scale = 0;
    NumericVector col_sum = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    if(left_index - 1 != 0){
      scale = col_sum[right_index - 1] - col_sum[left_index - 2];
      block_size = col_samples[right_index - 1] - col_samples[left_index - 2];
    } else {
      scale = col_sum[right_index - 1];
      block_size = col_samples[right_index - 1];
    }
    scale /= block_size;
    loglike = -block_size*(log(scale) + 1);
  }

  // log likelihood computation for Poisson distribution
  if(family == "poisson"){
    double rate = 0;
    NumericVector col_sum = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    if(left_index - 1 != 0){
      rate = col_sum[right_index - 1] - col_sum[left_index - 2];
      block_size = col_samples[right_index - 1] - col_samples[left_index - 2];
    } else {
      rate = col_sum[right_index - 1];
      block_size = col_samples[right_index - 1];
    }
    rate /= block_size;
    //not including factorial terms that do not depend on the rate
    loglike = block_size*rate*(log(rate) - 1);
  }


  return(-loglike);
}

float Blockcpd::compute_loss(const int& left_index,
                             const int& right_index){

  return compute_negloglike(left_index, right_index) + *REAL(pen_func(left_index, right_index));

}

float Blockcpd::compute_regularization(const int& left_index,
                                       const int& right_index){
 return  *REAL(pen_func(left_index, right_index));
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
    parameters = List::create(Named("prob") = probs);
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

  // Fitting parameters for binary Markov chain model
  if(family == "binaryMarkov"){
    NumericVector col_samples0 = suff_stats[0];
    NumericVector col_sums0 = suff_stats[1];
    NumericVector col_samples1 = suff_stats[2];
    NumericVector col_sums1 = suff_stats[3];
    int block_samples0;
    int block_samples1;
    NumericVector probs00(n_blocks);
    NumericVector probs01(n_blocks);
    NumericVector probs10(n_blocks);
    NumericVector probs11(n_blocks);
    for(int i = 0; i < n_blocks; i++){
      probs00[i] = 0;
      probs11[i] = 0;
      block_samples0 = 0;
      block_samples0 = 0;
      // Note 1 can not be a change point!
      if(block_begin[i] - 1 > 1){
        // the -3 on the left_index - 3 is due to the convention that
        // at the end of a block, it transits to the beginning of a new block
        // using the transition matrix of the next block!!!
        probs00[i] = col_sums0[block_end[i] - 2] - col_sums0[block_begin[i] - 3];
        probs11[i] = col_sums1[block_end[i] - 2] - col_sums1[block_begin[i] - 3];
        block_samples0 = col_samples0[block_end[i] - 2] - col_samples0[block_begin[i] - 3];
        block_samples1 = col_samples1[block_end[i] - 2] - col_samples1[block_begin[i] - 3];
      } else {
        probs00[i] = col_sums0[block_end[i] - 2];
        probs11[i] = col_sums1[block_end[i] - 2];
        block_samples0 = col_samples0[block_end[i] - 2];
        block_samples1 = col_samples1[block_end[i] - 2];
      }
      probs00[i] /= block_samples0;
      probs01[i] = 1 - probs00[i];
      probs11[i] /= block_samples1;
      probs10[i] = 1 - probs11[i];
    }
    parameters = List::create(Named("prob00") = probs00,
                              Named("prob01") = probs01,
                              Named("prob10") = probs10,
                              Named("prob11") = probs11);
  }
  // Fitting parameters for exponential distribution
  if(family == "exponential"){
    NumericVector col_sums = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    int block_size;
    NumericVector scales(n_blocks);
    for(int i = 0; i < n_blocks; i++){
      scales[i] = 0;
      block_size = 0;
      if(block_begin[i] - 1 != 0){
        scales[i] = col_sums[block_end[i] - 1] - col_sums[block_begin[i] - 2];
        block_size = col_samples[block_end[i] - 1] - col_samples[block_begin[i] - 2];
      } else {
        scales[i] = col_sums[block_end[i] - 1];
        block_size = col_samples[block_end[i] - 1];
      }
      scales[i] /= block_size;
    }
    parameters = List::create(Named("scale") = scales);
  }

  // Fitting parameters for Poisson distribution
  if(family == "poisson"){
    NumericVector col_sums = suff_stats[0];
    NumericVector col_samples = suff_stats[1];
    int block_size;
    NumericVector rates(n_blocks);
    for(int i = 0; i < n_blocks; i++){
      rates[i] = 0;
      block_size = 0;
      if(block_begin[i] - 1 != 0){
        rates[i] = col_sums[block_end[i] - 1] - col_sums[block_begin[i] - 2];
        block_size = col_samples[block_end[i] - 1] - col_samples[block_begin[i] - 2];
      } else {
        rates[i] = col_sums[block_end[i] - 1];
        block_size = col_samples[block_end[i] - 1];
      }
      rates[i] /= block_size;
    }
    parameters = List::create(Named("rate") = rates);
  }

}

void Blockcpd::sort_changepoints(){
  std::sort(changepoints.begin(), changepoints.end());
}

