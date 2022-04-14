#include <Rcpp.h>
#include "Hierseg.h"
#include <algorithm>
#include <queue>
using namespace Rcpp;

Hierseg::Hierseg(String family, const List& suff_stats, Function pen_func,
                 int ncol, int min_block_size, int max_blocks,
                 String algorithm_type)
  : Blockcpd(family, suff_stats, pen_func, ncol, min_block_size, max_blocks),
    algorithm_type(algorithm_type){}

// Wrap for the fitting process
void Hierseg::fit_hierseg(){
  //----
  // Fits the change point set
  float current_nll = compute_negloglike(1, ncol);
  float current_loss = compute_loss(1, ncol);
  if(algorithm_type == "recursive"){
    binary_split(1, ncol, current_nll, current_loss);
  }
  if(algorithm_type == "iterative"){
    binary_split_iter(current_nll,  current_loss);
  }
  sort_changepoints();
  //---
  // Fits family parameters given the change point set
  fit_family_parameters();
}

// Auxiliary function
bs_node Hierseg::get_best_split(const int& left_index,
                                const int& right_index){
  Rcpp::checkUserInterrupt(); // check user interruption in rcpp
  float initial_nll = compute_negloglike(left_index, right_index);
  float initial_loss = initial_nll + compute_regularization(left_index, right_index);
  float left_loss = 0;
  float right_loss = initial_loss;
  float left_nll = 0;
  float right_nll = initial_nll;
  int split_index = 0;
  float new_left_loss, new_right_loss;
  float new_left_nll, new_right_nll; //left and right negloglike
  bs_node node;
  node.left = left_index;
  node.right = right_index;
  //---
  // Forced halt conditions
  if(right_index - left_index + 1 < 2*min_block_size){
    node.split_index = 0;
    node.loss_reduction = 0;
    node.nll_reduction = 0;
    return node;
  } // if this restriction is true, we can not split further, so we halt early
  if(left_index == right_index) {
    node.split_index = 0;
    node.loss_reduction = 0;
    node.nll_reduction = 0;
    return node;
  } // this is a subcase of the above since min_block_size >= 1. Just to ensure.
  //---


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

  node.split_index = split_index;
  node.loss_reduction = initial_loss - (left_loss + right_loss);
  node.nll_reduction = initial_nll - (left_nll + right_nll);
  return node;
}

// Iterative Implementation
void Hierseg::binary_split_iter(const float& unsplit_nll,
                                const float& unsplit_loss) {
  // Final model neg loglike and loss
  // they are class attributes
  negloglike = unsplit_nll;
  loss = unsplit_loss;
  int n_changepoints = 0;

  // Setting up data structure
  std::priority_queue<bs_node> split_queue;
  bs_node search_node;
  bs_node curr_node;
  search_node = get_best_split(1, ncol);

  if(search_node.split_index != 0){
    split_queue.push(search_node);
    loss -= search_node.loss_reduction;
    negloglike -= search_node.nll_reduction;
  }

  // First interval
  while(!split_queue.empty()){
    // add split index to change points and remove from queue
    curr_node = split_queue.top();
    changepoints.push_back(curr_node.split_index);
    n_changepoints++;
    split_queue.pop();

    // halt condition
    if(n_changepoints >= max_blocks){
      // empties queue
      while(!split_queue.empty()){split_queue.pop();}
      return;
    }

    // two new intervals to search
    // search [left, split_index)
    search_node = get_best_split(curr_node.left,
                                 curr_node.split_index);
    if(search_node.split_index != 0){
      split_queue.push(search_node);
      loss -= search_node.loss_reduction;
      negloglike -= search_node.nll_reduction;
    }
    // search [split_index+1, right)
    search_node = get_best_split(curr_node.split_index + 1,
                                 curr_node.right);
    if(search_node.split_index != 0){
      split_queue.push(search_node);
      loss -= search_node.loss_reduction;
      negloglike -= search_node.nll_reduction;
    }
  }
  return;
}

// Recursive implementation
void Hierseg::binary_split(const int& left_index,
                            const int& right_index,
                            const float& current_nll,
                            const float& current_loss) {
  Rcpp::checkUserInterrupt(); // check user interruption in rcpp
  if(right_index - left_index + 1 < 2*min_block_size){
    loss += current_loss;
    negloglike += current_nll;
    return;
  } // if this restriction is true, we can not split further, so we halt early
  if(left_index == right_index) {
    loss += current_loss;
    negloglike += current_nll;
    return;
  } // this is a subcase of the above since min_block_size >= 1. Just to ensure.

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
