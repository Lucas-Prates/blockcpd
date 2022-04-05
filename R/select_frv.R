#' @title
#' Methodology to aid choosing regularization constant
#'
#' @description
#' Aids in the selection of the penalization constants, possibly providing an
#' automatic optimal value. It analyses how the number of change-points vary
#' with the chosen grid of penalization constant. It applied the
#' First Repeated Value (FRV) methodology
#' to select the regularization constant lambda. It is similar to the Elbow
#' method used in clustering, or the CROPS algorithm in change-point detection.
#' The values of the constant range from 'lambda_left' to
#' 'lambda_right', increasing by 'step'. For each value, the function
#' \link[=fit_blockcpd]{fit_blockcpd} is run with arguments 'model_args'. An
#' automatic suggestion for the penalization, number of change-points and model
#' is given automatically. Optionally, The user can call the plot function to
#' the output of this method so he can use an elbow plot like graphical
#' inspection to select the constant value.
#'
#' @param data_matrix Data frame or matrix containing the data set to be
#' segmented.
#' @param lambda_left Left most value of lambda. Must be non-negative.
#' @param lambda_right Right most value of lambda. Must be non-negative and
#' greater than lambda_left.
#' @param step Value by which lambda will be increased. Must be greater than 0,
#' The default is 'automatic', which consists of a penalization of
#' 1/sqrt(log(n)), where n is the number of samples (rows).
#' @param model_args A list with argument values for the
#' \link[=fit_blockcpd]{fit_blockcpd} function. The list keys must be the
#' arguments names. It must *not* contain the argument 'lambda' or
#' 'data_matrix'.
#'
#' @return Returns a frv object containing the suggested values and caller
#' parameters.
#' @export
select_frv = function(data_matrix, lambda_left = 0, lambda_right = 10,
                      step = "automatic", model_args = list()){

  # check input
  if(step == "automatic"){
    n = nrow(data_matrix)
    m = ncol(data_matrix)
    step = 1/sqrt(log(n*m))
  }
  args_to_check = list(model_args = model_args,
                       lambda_left = lambda_left,
                       lambda_right = lambda_right,
                       step = step)
  do.call(check_input, list(caller = "select_frv",
                            args_to_check = args_to_check))


  lambda_set = seq(lambda_left, lambda_right, step)
  lambda_set_len = length(lambda_set)
  # guarantees that there are at least 5 values in the grid
  if(lambda_set_len < 5){
    warning("Fewer than 5 constant values in the specified range. Changing the
            range to have at least 5 points!")
    lambda_set = seq(lambda_left, lambda_right, length = 5)
    lambda_set_len = 5
  }
  suggested_lambda = NULL
  suggested_ncp = NULL
  ncp = numeric(lambda_set_len)
  neg_loglike = numeric(lambda_set_len)
  call_arg = c(list(data_matrix = data_matrix, lambda = NULL),
               model_args)
  call_arg$skip_input_check = TRUE # skips input to avoid error
  previous_ncp = ncol(data_matrix) + 1 # impossible
  for(i in 1:lambda_set_len){
    call_arg$lambda = lambda_set[i]
    model = do.call(fit_blockcpd, call_arg)
    ncp[i] = model$ncp
    neg_loglike[i] = model$neg_loglike
    # FRV: choose the first value in which the estimated ncp repeats its value
    if((ncp[i] == previous_ncp)&&(is.null(suggested_lambda))){
      suggested_lambda = lambda_set[i]
      suggested_ncp = ncp[i]
      suggested_model = model
    }
    previous_ncp = ncp[i]
  }
  # If no repetition appeared, we estimate the best constant by estimating the
  # curvature. This is not the FRV per se, but its a fast workaround.
  # TO BE DONE: implement the FRV properly, refining the grid instead of using
  # curvature estimation
  if(is.null(suggested_lambda)){
    curvature_info = compute_max_curvature(lambda_set, ncp, step)
    suggested_lambda = curvature_info$argmax_curv
    suggested_ncp = curvature_info$fmax_curv
    call_arg$lambda = suggested_lambda
    suggested_model = do.call(fit_blockcpd, call_arg)
  }
  call_data = list(lambda_left = lambda_left, lambda_right = lambda_right,
                   step = step, model_args = model_args)
  frv = list(lambda_set = lambda_set, ncp = ncp,
                  neg_loglike = neg_loglike,
                  suggested_lambda = suggested_lambda,
                  suggested_ncp = suggested_ncp,
                  suggested_model = suggested_model,
                  call_data = call_data
                  )
  class(frv) = "frv"
  return(frv)
}
