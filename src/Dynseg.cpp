#include <Rcpp.h>
#include "Dynseg.h"
#include <unistd.h>
using namespace Rcpp;

Dynseg::Dynseg(String family, const List& suff_stats, Function pen_func,
               int ncol, int min_block_size, int max_blocks)
  : Blockcpd(family, suff_stats, pen_func, ncol, min_block_size),
    max_blocks(max_blocks),
    loss_mat(ncol, ncol) {}

void Dynseg::fit_dynseg(){
  build_loss_mat();
  fit_cp_set();
  fit_family_parameters();

  // Computes neg log-likelihood
  negloglike = loss;
  if(changepoints.size() > 0){
    negloglike -= *REAL(pen_func(1, changepoints[0]));
    for(int i = 0; i < changepoints.size() - 1; i++){
      negloglike -= *REAL(pen_func(changepoints[i] + 1, changepoints[i + 1]));
    }
    negloglike -= *REAL(pen_func(changepoints[changepoints.size()] + 1, ncol));
  }
  else{
    negloglike -= *REAL(pen_func(1, ncol));
  }

}

// !!! Lazy implementation !!!
// A more efficient algorithm can be written to compute the loss,
// using a recursion that depends on the family.
void Dynseg::build_loss_mat(){

  for(int i = 0; i < ncol; i++){
    Rcpp::checkUserInterrupt();
    for(int j = i; j < ncol; j++){
      // instead of calling the function, an efficient approach would
      // compute the loss using dynamical programming.
      loss_mat.set_value(compute_loss(i + 1, j + 1), i, j);
    }
  }
}


// Pseudocode provided at paper **put link to paper**
void Dynseg::fit_cp_set(){
  // The element best_loss[k][j] refers to the best loss from 1 to j with
  // exactly k change points.
  TriangularMatrix<double> best_loss(max_blocks + 1, ncol);

  // The element cp_mat[k][j] refers to the location of the k-th change point
  // considering the best split of the interval [1, j] with exactly k change
  // points
  TriangularMatrix<int> cp_mat(max_blocks + 1, ncol);

  for(int j = 0; j < ncol; j++){
    //best_loss[0][j] = loss_mat[0][j]
    best_loss.set_value(loss_mat.get_value(0, j), 0, j);
    cp_mat.set_value(ncol - 1, 0, j);
  }
  int best_k = 0; // best number of change points
  float min_loss = best_loss.get_value(0, ncol - 1); // loss for best_k

  int iter_best_cp; // best location for cp split at each iteration
  double iter_min_loss; // best loss value for each iteration
  double curr_loss; // loss at current step

  if (max_blocks >= 1){
    for (int k = 1; k <= max_blocks; k++){
      checkUserInterrupt();
      for (int j = k; j < ncol; j++){
        iter_min_loss = best_loss.get_value(k - 1, k - 1) + loss_mat.get_value(k, j);
        iter_best_cp = k;
        for (int i = k; i <= j; i++){
          curr_loss = best_loss.get_value(k - 1, i - 1) + loss_mat.get_value(i, j);
          // update local (iteration) loss and change point number
          if(iter_min_loss > curr_loss){
            iter_min_loss = curr_loss;
            iter_best_cp = i;
          }
        }
        best_loss.set_value(iter_min_loss, k, j);
        cp_mat.set_value(iter_best_cp, k, j);
      }
      // update global loss and change point number
      if (min_loss > best_loss.get_value(k, ncol - 1)){
        min_loss = best_loss.get_value(k, ncol - 1);
        best_k = k;
      }
    }
  }
  loss = min_loss; // update class loss
  // Retrieve change point set
  int segf, segi;

  segf = ncol - 1;
  segi = 0;

  if (best_k > 0) {
    for (int j = best_k; j > 0; j--){
      segi = cp_mat.get_value(j, segf);
      changepoints.insert(changepoints.begin(), segi);
      segf = segi - 1;
    }
  }
}
