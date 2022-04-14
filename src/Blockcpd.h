#ifndef BLOCKCPD_H
#define BLOCKCPD_H

#include <Rcpp.h>
using namespace Rcpp;

// Class for the statistical model used to fit the data. It contains the
// family, sufficient statistics, regularization function, loss, negative
// log-likelihood and parameters. It provides methods to compute the negloglike,
// regularization loss and fit family parameters. However, it does not implement
// a method to estimate the change point set. This is implemented in the
// subclasses Hierseg and Dynseg.
class Blockcpd{
public:
  String family; // family fitted
  List suff_stats; // list of sufficient statistics for the family
  Function pen_func; // penalization function defined by user
  int ncol;
  int min_block_size;
  int max_blocks;
  std::vector<int> changepoints;
  float loss; // total regularized loss function for the estimated model
  float negloglike; // total negative log-likelihood
  List parameters;

  Blockcpd(String family, const List& suff_stats,
           Function pen_func, int ncol,
           int min_block_size, int max_blocks);

  // Computes the negative log-likelihood for the block defined by the indices
  // Indices are passed considering start index as 1 (not 0, as usual in c++)
  float compute_negloglike(const int& left_index, const int& right_index);

  // Returns the regularized loss of the model. It is given by the negloglike
  // plus the penalization function provided by the user
  float compute_loss(const int& left_index, const int& right_index);

  // Returns regularization value on interval
  float compute_regularization(const int& left_index, const int& right_index);

  void sort_changepoints();

  // FITS -> Parameters
  // fit family parameters AFTER the change point set has been fitted
  void fit_family_parameters();

};

#endif
