#' @title Implements the regularization functions used in the estimation
#'
#' @description The estimator in this package computes the optimum of
#' \eqn{-l(C, p) + \lambda*R(leftIndex, rightIndex, n, m)}, where l is the
#' log likelihood of the family, lambda is the penalization constant and R
#' is the regularization function. The user can create his own regularization
#' function and pass as an argument to \link[fit_blockcpd]{fit_blockcpd}. It
#' should have four arguments, in the following order: leftIndex, rightIndex, n
#' (the sample size) and m (the number of columns or variables).
#' If the function depends on leftIndex and rightIndex, it will be
#' non-homogeneous, which might be interesting in some applications.
#' The package implements some functions as an example, but uses only
#' bic_loss_hs as the default. The algorithm is consistent as long as the
#' the regularization is bounded by a constant.
#'
#' @examples
#' my_reg <- function(leftIndex, rightIndex, n, m){
#'     block_size = (rightIndex - leftIndex + 1)
#'     return(log(n*m)*(1/block_size))
#' }
#' @rdname regularization
toy_regularization = function(left_index, right_index, n, m){
  block_size = (rightIndex - leftIndex + 1)
  return(log(n*m)*(1/block_size))
}

# Regular BIC loss
bic_loss = function(left_index, right_index, n, m){
  loss = log(n)
  return(loss)
}

# Square root loss in n
sqrt_loss = function(left_index, right_index, n, m){
  loss = sqrt(n)
  return(loss)
}

# Log loss on n with exponential decay on block size. This favors bigger
# blocks.
exp_loss = function(left_index, right_index, n, m){
  loss = log(n)*(1+exp(-(right_index - left_index)))
  return(loss)
}

sqrt_n_exp_loss = function(left_index, right_index, n, m){
  loss = sqrt(n)*(1+exp(-(right_index - left_index)))
  return(loss)
}
