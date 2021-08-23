#' @title
#' Block segmentation using hierarchical algorithm
#'
#' @description
#' Uses binary splitting to obtain a greedy soltuion to the regularized loss
#' optimization problem. Should be called within
#' \link[=fit_blockcpd]{fit_blockcpd}
#'
#' @param suff_stats Sufficient statistics to perform change point analysis
#' @param family The name of the family used to fit the model
#' @param lambda Penalization constant
#' @param n Number of samples
#' @param m Number of variables
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, n, m),
#'   where the left_index:right_index is the integer interval, n the sample
#'   size and m the number of variables/columns.
compute_hierseg = function(suff_stats,
                           family,
                           lambda = 1,
                           n,
                           m,
                           pen_func = bic_loss,
                           max_blocks = NULL) {

  # max_blocks is not used in this function

  # Penalization function that will be called in .cpp extension
  hs_pen_function = function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, n, m) )
  }

  hs_output = compute_hierseg_cpp(suff_stats = suff_stats,
                                  family = family,
                                  ncol = m,
                                  hs_pen_function)

  model_info = list(changepoints = hs_output[[1]],
                    parameters = hs_output[[2]],
                    loss = hs_output[[3]],
                    neg_loglike = hs_output[[4]],
                    n_cp = length(hs_output[[1]])
                    )

  return(model_info)

}
