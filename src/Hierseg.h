#ifndef HIERSEG_H
#define HIERSEG_H

#include <Rcpp.h>
#include "Blockcpd.h"
using namespace Rcpp;

// Auxiliar data structure for fitting binary split iteratively
typedef struct{
  unsigned int split_index;
  unsigned int left, right;
  float nll_gain; // negative log likelihood
  float loss_gain;
} bs_node;

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

  bs_node* get_best_split(unsigned int left_index,
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
  // Perform calculations, selection best splitting indexes and appending to
  // the change point set of the stats model
  void binary_split_iter(const float& unsplit_nll,
                         const float& unsplit_loss);
};

#endif
