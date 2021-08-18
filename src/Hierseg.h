#ifndef HIERSEG_H
#define HIERSEG_H

#include <Rcpp.h>
#include "Blockcpd.h"
using namespace Rcpp;

// Class used to fit the data using the hierarchical algorithm. Inherits from
// Blockcpd, which provides the members and methods of the statistical model.
// The importance of this class is to estimate the change point set.
class Hierseg : public Blockcpd
{
public:

  Hierseg(String family, const List& suff_stats, Function pen_func, int ncol);

  // Wrapper for fitting methods. First, it calls a method to fit the change
  // point set. Then, it calls fit_family_parameters.
  void fit_hierseg();


  // FITS -> Change point set
  // Recursive implementation of the hierarchical algorithm
  // Perform calculations, selection best splitting indexes and appending to
  // the change point set of the stats model
  void binary_split(const int& left_index,
                    const int& right_index,
                    const float& current_nll,
                    const float& current_loss);
};

#endif
