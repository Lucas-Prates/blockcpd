#ifndef TRIANGULAR_MATRIX_H
#define TRIANGULAR_MATRIX_H

#include <Rcpp.h>

// Triangular matrix struct used in the dynamical programming algorithm
template <typename T>
class TriangularMatrix{
public:
  int nrow, ncol;
  int total_size;
  std::vector<T> tmat;

  TriangularMatrix(int nrow, int ncol);

  double get_value(int i, int j);

  void set_value(double value, int i, int j);

};

#endif
