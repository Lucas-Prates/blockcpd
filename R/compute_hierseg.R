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
#' @param nrow Number of rows or samples
#' @param ncol Number of columns or variables
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, nrow, ncol),
#'   where the left_index:right_index is the integer interval, nrow the sample
#'   size and ncol the number of variables/columns.
#' @param min_block_size Minimum block size allowed. Default is 0, and the value
#' must be smaller or equal to ncol.
#' @param max_blocks Threshold on the number of block segments to fit the model.
#' Set low values for this parameters if having performance issues on large
#' data sets.
compute_hierseg = function(suff_stats,
                           family,
                           lambda = 1,
                           nrow,
                           ncol,
                           pen_func = bic_loss,
                           min_block_size = min_block_size,
                           max_blocks = NULL) {

  if(is.null(max_blocks)){
    max_blocks = ncol-1
  }

  # Penalization function that will be called in .cpp extension
  hs_pen_function = function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, nrow, ncol) )
  }

  hs_output = compute_hierseg_cpp(suff_stats = suff_stats,
                                  family = family,
                                  ncol = ncol,
                                  min_block_size = min_block_size,
                                  max_blocks = max_blocks,
                                  pen_func = hs_pen_function,
                                  algorithm_type = "iterative")

  model_info = list(changepoints = hs_output[[1]],
                    parameters = hs_output[[2]],
                    loss = hs_output[[3]],
                    neg_loglike = hs_output[[4]],
                    ncp = length(hs_output[[1]])
                    )

  return(model_info)

}
