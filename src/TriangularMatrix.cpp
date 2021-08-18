#include <Rcpp.h>
using namespace Rcpp;

// Triangular matrix struct used in the dynamical programming algorithm
template <typename T>
class TriangularMatrix{
public:
  int nrow, ncol;
  int total_size;
  std::vector<T> tmat;

  TriangularMatrix(int nrow, int ncol)
    : nrow(nrow), ncol(ncol)
  {
    if(nrow > ncol){
      throw std::invalid_argument("TriangularMatrix must have nrow <= ncol, received nrow > ncol");
    }
    total_size = (ncol*(ncol + 1))/2 - (((ncol - nrow)*(ncol - nrow + 1))/2);
    std::vector<T> aux_T(total_size);
    tmat = aux_T;
  }

  double get_value(int i, int j){
    if(i > nrow){
      throw std::invalid_argument("TriangularMatrix: i > nrow, acessing non-existing row");
    }
    if(j < i){
      throw std::invalid_argument("TriangularMatrix: j must be greater or equal to i");
    }
    int offset = (ncol + 1)*i - (i*(i + 1))/2;
    if(offset + j - i > total_size){
      printf("TriangularMatrix: index_sum: %d, total_size: %d\n", offset + j - i, total_size);
      throw std::invalid_argument("TriangularMatrix: accessing value out of range");
    }
    return(tmat[offset + j - i]);
  }

  void set_value(double value, int i, int j){
    if(i > nrow){
      throw std::invalid_argument("TriangularMatrix: i > nrow, acessing non-existing row");
    }
    if(j < i){
      throw std::invalid_argument("TriangularMatrix: j must be greater or equal to i");
    }
    int offset = (ncol + 1)*i - (i*(i + 1))/2;
    if(offset + j - i > total_size){
      printf("TriangularMatrix: index_sum: %d, total_size: %d\n", offset + j - i, total_size);
      throw std::invalid_argument("TriangularMatrix: accessing value out of range");
    }
    tmat[offset + j - i] = value;
  }

  };
