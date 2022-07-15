#ifndef DYNSEG_H
#define DYNSEG_H

#include <Rcpp.h>
#include "TriangularMatrix.h"
#include "Blockcpd.h"
using namespace Rcpp;

// Class used to fit the data using the dynamical programming
// algorithm. Inherits from Blockcpd, which provides the members and methods of
// the statistical model.
// The importance of this class is to estimate the change point set.
class Dynseg : public Blockcpd{
public:

  // The element loss_mat[i][j] corresponds to the loss of the interval [i, j]
  TriangularMatrix<double> loss_mat;

  Dynseg(String family, const List& suff_stats, Function pen_func,
         int ncol, int min_block_size, int max_blocks);

  // Wrapper for fitting methods. First, it calls a method to fit the change
  // point set. Then, it calls fit_family_parameters.
  void fit_dynseg();

  // Build DP recursion matrix
  void build_loss_mat();

  // Fits the cp set using DP
  void fit_cp_set();

};

#endif
