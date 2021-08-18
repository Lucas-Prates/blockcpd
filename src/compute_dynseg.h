#ifndef COMPUTE_DYNSEG_H
#define COMPUTE_DYNSEG_H

#include <Rcpp.h>
using namespace Rcpp;


List compute_dynseg_cpp(const NumericMatrix& data_mat,
                        const int& ncol,
                        int segthr,
                        const Function& pen_func);


#endif
