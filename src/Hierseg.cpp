#include <Rcpp.h>
#include "Hierseg.h"
using namespace Rcpp;

Hierseg::Hierseg(String family, const List& suff_stats, Function pen_func,
                 int ncol, int min_block_size)
  : Blockcpd(family, suff_stats, pen_func, ncol, min_block_size) {}


void Hierseg::fit_hierseg(){
  //----
  // Fits the change point set
  float current_nll = compute_negloglike(1, ncol);
  float current_loss = compute_loss(1, ncol);
  binary_split(1, ncol, current_nll, current_loss);

  sort_changepoints();
  //---

  // Fits family parameters given the change point set
  fit_family_parameters();

}

void Hierseg::binary_split(const int& left_index,
                            const int& right_index,
                            const float& current_nll,
                            const float& current_loss) {

  Rcpp::checkUserInterrupt(); // check user interruption in rcpp
  if(left_index == right_index) {
    loss += current_loss;
    negloglike += current_nll;
    return;
  }

  float left_loss = 0;
  float right_loss = current_loss;
  float left_nll = 0;
  float right_nll = current_nll;
  int split_index = 0;
  float new_left_loss, new_right_loss;
  float new_left_nll, new_right_nll; //left and right negloglike

  for(int j = left_index; j < right_index; j++) {

    new_left_nll = compute_negloglike(left_index, j);
    //negloglike + regularization
    new_left_loss = new_left_nll + compute_regularization(left_index, j);

    new_right_nll = compute_negloglike(j + 1, right_index);
    //negloglike + regularization
    new_right_loss = new_right_nll + compute_regularization(j + 1, right_index);

    if(new_left_loss + new_right_loss < left_loss + right_loss){
      left_loss = new_left_loss;
      right_loss = new_right_loss;
      left_nll = new_left_nll;
      right_nll = new_right_nll;
      split_index = j;
    }


  }

  if(split_index != 0) {
    changepoints.push_back(split_index);
    binary_split(left_index, split_index, left_nll, left_loss);

    binary_split(split_index + 1, right_index, right_nll, right_loss);
  } else {
    loss += current_loss;
    negloglike += current_nll;
  }

  return;

}

