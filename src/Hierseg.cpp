#include <Rcpp.h>
#include "Hierseg.h"
#include <algorithm>
#include <utility>
using namespace Rcpp;

Hierseg::Hierseg(String family, const List& suff_stats, Function pen_func,
                 int ncol, int min_block_size)
  : Blockcpd(family, suff_stats, pen_func, ncol, min_block_size) {}


void Hierseg::fit_hierseg(){
  //----
  // Fits the change point set
  float current_nll = compute_negloglike(1, ncol);
  float current_loss = compute_loss(1, ncol);
  //binary_split(1, ncol, current_nll, current_loss);
  binary_split_iter(current_nll,  current_loss);
  sort_changepoints();
  //---
  // Fits family parameters given the change point set
  fit_family_parameters();

}

// Auxiliary function for fitting methods
// ---
// Input
// left_index: (unsigned int) left index of the interval
// right_index: (unsigned int) left index of the interval
// ---
// Output
// node: (bs_node*) Returns a pointer to a structure that provides
// the best split point, the left and right indices from the call,
// and the log gains of the negloglike and loss
bs_node* Hierseg::get_best_split(const unsigned int& left_index,
                                const unsigned int& right_index){
  Rcpp::checkUserInterrupt(); // check user interruption in rcpp
  float initial_loss = compute_loss(left_index, right_index);
  float initial_nll = compute_negloglike(left_index, right_index);
  float left_loss = 0;
  float right_loss = initial_loss;
  float left_nll = 0;
  float right_nll = initial_nll;
  int split_index = 0;
  float new_left_loss, new_right_loss;
  float new_left_nll, new_right_nll; //left and right negloglike
  bs_node* node;
  node->left = left_index;
  node->right = right_index;
  //---
  // Forced halt conditions
  if(right_index - left_index + 1 < 2*min_block_size){
    node->split_index = 0;
    node->loss_gain = 0;
    node->nll_gain = 0;
    return node;
  } // if this restriction is true, we can not split further, so we halt early
  if(left_index == right_index) {
    node->split_index = 0;
    node->loss_gain = 0;
    node->nll_gain = 0;
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

  node->split_index = split_index;
  node->loss_gain = initial_loss - (left_loss + right_loss);
  node->nll_gain = initial_nll - (left_nll + right_nll);
  return node;
}

// Iterative Implementation
// The node_gain pair consists of a the node split info + the log gain
// the nodes are sorted by their neg loss gain value. Nodes with greater neg gain
// are added first;
typedef std::pair<bs_node*, float> node_gain;
void Hierseg::binary_split_iter(const float& unsplit_nll,
                                const float& unsplit_loss) {
  // Final model neg loglike and loss
  // they are class attributes
  negloglike = unsplit_nll;
  loss = unsplit_loss;

  // Setting up data structure
  std::vector<node_gain> split_queue;
  std::make_heap(split_queue.begin(), split_queue.end());
  bs_node* search_node;
  node_gain curr_node;
  search_node = get_best_split(1, ncol);
  float loss_gain, nll_gain;
  if(search_node->split_index != 0){
    split_queue.push_back(node_gain(search_node, search_node->loss_gain));
    std::push_heap(split_queue.begin(), split_queue.end());
    loss += search_node->loss_gain;
    negloglike += search_node->nll_gain;
  }

  // First interval
  Rcout << "Reached iterations for split_queue\n";
  while(!split_queue.empty()){
    // add split index to change points and remove from queue
    curr_node = split_queue.front();
    changepoints.push_back(curr_node.first->split_index);
    std::pop_heap(split_queue.begin(), split_queue.end());
    split_queue.pop_back();
    Rcout << "Added " << curr_node.first->split_index << " to changepoints\n";

    // two new intervals to search
    // search [left, split_index)
    search_node = get_best_split(curr_node.first->left,
                                 curr_node.first->split_index);
    if(search_node->split_index != 0){
      split_queue.push_back(node_gain(search_node, search_node->loss_gain));
      std::push_heap(split_queue.begin(), split_queue.end());
      loss += search_node->loss_gain;
      negloglike += search_node->nll_gain;
    }
    // search [split_index+1, right)
    search_node = get_best_split(curr_node.first->split_index + 1,
                                 curr_node.first->right);
    if(search_node->split_index != 0){
      split_queue.push_back(node_gain(search_node, search_node->loss_gain));
      std::push_heap(split_queue.begin(), split_queue.end());
      loss += search_node->loss_gain;
      negloglike += search_node->nll_gain;
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
