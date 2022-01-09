#ifndef HIERSEG_H
#define HIERSEG_H

#include <Rcpp.h>
#include "Blockcpd.h"
using namespace Rcpp;

//---
// Auxiliary data structure for fitting binary split iteratively
// It contains the information that goes to the priority queue.
class bs_node{
public:
  unsigned int split_index;
  unsigned int left, right;
  float nll_reduction; // negative log likelihood
  float loss_reduction;
};


// overload for '<' operator to correctly compare nodes
inline bool operator<(const bs_node& node1, const bs_node& node2){
 return node1.loss_reduction < node2.loss_reduction;
}
//---

// Class used to fit the data using the hierarchical algorithm. Inherits from
// Blockcpd, which provides the members and methods of the statistical model.
// The importance of this class is to estimate the change point set.
class Hierseg : public Blockcpd
{
private:
  String algorithm_type; // A string to decide to use recursive or iterative
                         // implementation of the binary split algorithm
public:
  Hierseg(String family, const List& suff_stats, Function pen_func,
          int ncol, int min_block_size, int max_blocks, String algorithm_type);

  // Wrapper for fitting methods. First, it calls a method to fit the change
  // point set. Then, it calls fit_family_parameters.
  void fit_hierseg();

  // Auxiliary function for fitting methods
  // Given the end points of the search interval, it returns a bs_node object
  // containing the required information for fitting the binary segmentation
  // iteratively
  bs_node get_best_split(unsigned int left_index,
                          unsigned int right_index);

  // FITS -> Change point set
  // Recursive implementation of the hierarchical algorithm
  // Perform calculations, selection best splitting indexes and appending to
  // the change point set of the stats model
  void binary_split(const int& left_index,
                    const int& right_index,
                    const float& current_nll,
                    const float& current_loss);

  // FITS -> Change point set
  // Iterative implementation of the hierarchical algorithm
  // Uses a priority queue of bs_node so that change-points are added by a
  // measure of the loss reduction. This allows objective usage of the
  // max_blocks argument, as it guarantees that the added change-points are the
  // "best"
  void binary_split_iter(const float& unsplit_nll,
                         const float& unsplit_loss);
};

#endif
