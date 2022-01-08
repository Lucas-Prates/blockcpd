#' @title
#' Block segmentation using dynamical programming
#'
#' @description
#' Computes the exact solution of the regularized loss optimization problem,
#' providing change point locations and the parameters of each blocks. Should
#' be called within \link[=fit_blockcpd]{fit_blockcpd}
#'
#' @param suff_stats Sufficient statistics to perform change point analysis
#' @param family The name of the family used to fit the model
#' @param lambda Penalization constant
#' @param nrow Number of rows or samples
#' @param ncol Number of columns or variables
#' @param min_block_size Minimum block size allowed. Default is 0, and the value
#' must be smaller or equal to ncol.
#' @param max_blocks Threshold on the number of block segments to fit the model.
#' Set low values for this parameters if having performance issues on large
#' data sets.
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, nrow, ncol),
#'   where the left_index:right_index is the integer interval, nrow the sample
#'   size and ncol the number of variables/columns.
compute_dynseg = function(suff_stats,
                          family,
                          lambda = 1,
                          nrow,
                          ncol,
                          min_block_size = min_block_size,
                          max_blocks = ncol - 1,
                          pen_func = bic_loss){

  if (is.null(max_blocks)) {max_blocks = ncol - 1}
  if (max_blocks >= ncol) {max_blocks = ncol - 1}
  if (family == "binaryMarkov"){max_blocks = min(floor(ncol/2), max_blocks)}
  if (max_blocks <= 0) {max_blocks = 1}

  # Penalization function that will be called in .cpp extension
  ds_pen_function = function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, nrow, ncol) )
  }

  ds_output = compute_dynseg_cpp(suff_stats = suff_stats,
                                 family = family,
                                 ncol = ncol,
                                 min_block_size = min_block_size,
                                 max_blocks = max_blocks,
                                 pen_func = ds_pen_function)

  model_info = list(changepoints = ds_output[[1]],
                    parameters = ds_output[[2]],
                    loss = ds_output[[3]],
                    ncp = ds_output[[4]],
                    neg_loglike = ds_output[[5]]
                    )

  return(model_info)

}
