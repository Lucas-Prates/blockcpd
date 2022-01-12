#' @title
#' Plot for blockcpd object
#'
#' @description
#' Plots the selected parameters in a blocked fashion.
#'
#' @param blockcpd_obj A fitted blockcpd S3 object provided by the
#' \link[=fit_blockcpd]{fit_blockcpd} function.
#' @param parameter The parameter of the family for which to plot the blocked
#' @param index_values A numerical vector of size ncol that contains the values
#' of the the variable corresponding to the change points. For example, if your
#' segmented variable corresponds to a time samples from 0 to 150 sampled each
#' 15 seconds, the model treats these as values from 1 to 11. To plot on the
#' variable scale, pass the argument 'index_values = seq(0, 150, 15)'.
#' @param index_variable_name Name of the variable segmented.
#' @param pkg Graphical package to be used for plotting. Current values are
#' "base".
#'
#' @export
plot.blockcpd = function(blockcpd_obj, parameter = NULL,
                         index_values = NULL, index_variable_name = "Index",
                         pkg = "base"){
  # check if parameter argument is in family parameter list
  if(is.null(parameter)){
    parameter = names(blockcpd_obj$parameters)[1]
  }
  ncol = blockcpd_obj$metadata$ncol

  if(is.null(index_values)){
    index_values = 1:ncol
  }

  args_to_check = list(parameter = parameter,
                       family_parameters = names(blockcpd_obj$parameters),
                       is_index_values_numeric = is.numeric(index_values),
                       length_index_values = length(index_values),
                       ncol = ncol)
  do.call(check_input, list(caller = as.character(match.call()[[1]]),
                            args_to_check = args_to_check))

  ncp = blockcpd_obj$ncp
  changepoints = blockcpd_obj$changepoints
  parameter_vec = blockcpd_obj$parameters[[parameter]]
  if(pkg == "base"){
    sf = stepfun(index_values[changepoints], parameter_vec)
    plot(sf, xlim = c(index_values[1], index_values[ncol]), do.points = F,
         xaxs = "i", verticals = F,
         xlab = index_variable_name, ylab = parameter,
         main = paste("Block plot for", parameter, "parameter"))
    abline(v = blockcpd_obj$changepoints, col = "red", lty = "dashed")
  }

}

#' @title
#' Plot to aid choosing regularization constant
#'
#' @description
#' Plots how the number of change-points estimated by the given model vary
#' with the regularization constant lambda. It is similar to the Elbow method
#' used in clustering. The values of the constant range from 'lambda_left' to
#' 'lambda_right', increasing by 'step'. For each value, the function
#' \link[=fit_blockcpd]{fit_blockcpd} is run with arguments 'model_args'. The
#' modeler should choose the first value by which the number of change-points
#' curve starts to "flat-out" after that value. If there are at least 10 values
#' of constants to evaluate, a suggestion for lambda is also provided.
#'
#' @param data_matrix Data frame or matrix containing the data set to be
#' segmented.
#' @param lambda_left Left most value of lambda. Must be non-negative.
#' @param lambda_right Right most value of lambda. Must be non-negative and
#' greater than lambda_left.
#' @param step Value by which lambda will be increased. Must be greater than 0,
#' with default 1.
#' @param model_args A list with argument values for the
#' \link[=fit_blockcpd]{fit_blockcpd} function. The list keys must be the
#' arguments names. It must *not* contain the argument 'lambda' or
#' 'data_matrix'.
#' @param pkg Graphical package to be used for plotting. Current values are
#' "base".
#'
#' @return Along with the plot, it returns a list containing the lambda values,
#' number of change points per lambda, the negative log likelihood per lambda
#' and the model_args.
#' @export
elbow_plot = function(data_matrix, lambda_left = 0, lambda_right = 10,
                      step = 0.5, model_args = list(), pkg = "base"){

  # check input
  args_to_check = list(model_args = model_args,
                       lambda_left = lambda_left,
                       lambda_right = lambda_right,
                       step = step)
  do.call(check_input, list(caller = as.character(match.call()[[1]]),
                            args_to_check = args_to_check))


  lambda_set = seq(lambda_left, lambda_right, step)
  lambda_set_len = length(lambda_set)
  suggested_lambda = NULL
  suggested_ncp = NULL
  ncp = numeric(lambda_set_len)
  neg_loglike = numeric(lambda_set_len)
  call_arg = c(list(data_matrix = data_matrix, lambda = lambda_set[1]),
               model_args)
  call_arg$skip_input_check = TRUE # skips input to avoid error
  for(i in 1:lambda_set_len){
    call_arg$lambda = lambda_set[i]
    model = do.call(fit_blockcpd, call_arg)
    ncp[i] = model$ncp
    neg_loglike[i] = model$neg_loglike
  }
  # suggests a lambda based on estimated curvature
  if(lambda_set_len >= 10){
    # If there is a ncp value that repeats for different lambdas, choose the
    # smallest lambda such that it occurs.
    rep_index = which(diff(ncp) == 0)[1]
    if(!is.na(rep_index)){
      suggested_lambda = lambda_set[rep_index]
      suggested_ncp = ncp[rep_index]
    }
    # otherwise, make a suggestion based on estimated ncp curvature
    else{
      curvature_info = compute_max_curvature(lambda_set, ncp, step)
      suggested_lambda = curvature_info$argmax_curv
      suggested_ncp = curvature_info$fmax_curv
    }
  }
  if(pkg == "base"){
    plot(lambda_set, ncp,
         xlab = "lambda", ylab = "Number of change points",
         main = "Number of change points vs penalization constant")
    lines(lambda_set, ncp)
    if(!is.null(suggested_lambda)){
      abline(v = suggested_lambda, col = 'red', lty = "dashed")
      abline(h = suggested_ncp, col = 'red', lty = "dashed")
    }
  }
  elbow_plot_info = list(lambda = lambda_set, ncp = ncp,
                       neg_loglike = neg_loglike,
                       suggested_lambda = suggested_lambda,
                       suggested_ncp = suggested_ncp)

  invisible(elbow_plot_info)
}

#' @title
#' Plot to check reported change-points
#'
#' @description
#' Plots   the estimates of how likely it is for the model to detect a change at
#' any given point. True change-points should have confidence near $100%$,
#' while non change-points should have a confidence near $0%$. It might also be
#' difficult to detect a true change-point at the given sample size. In this
#' case, it should fluctuate in the middle.
#' @param model A `blockcpd` model object.
#' @param scale A string describing the scale which the y-scale should is
#' plotted. Possible values are "percentage", "probability" and "frequency".
#' @param index_values A numerical vector of size ncol that contains the values
#' of the the variable corresponding to the change points.
#' @param index_variable_name Name of the variable segmented.
#' @param pkg Graphical package to be used for plotting. Current values are
#' "base".
#'
#' @export
confidence_plot = function(model, scale = "percentage",
                           index_values = NULL, index_variable_name = "Index",
                           pkg = "base"){

  if(is.null(index_values)){index_values = 1:model$metadata$ncol}

  # check input
  args_to_check = list(model = model,
                       scale = scale,
                       is_index_values_numeric = is.numeric(index_values),
                       length_index_values = length(index_values),
                       ncol = model$metadata$ncol)
  do.call(check_input, list(caller = as.character(match.call()[[1]]),
                            args_to_check = args_to_check))

  scale_multiplier = 1
  n_samples = model$bootstrap_info$bootstrap_samples
  if(scale == "percentage"){scale_multiplier = 100/n_samples }
  if(scale == "probability"){scale_multiplier = 1/n_samples}
  cp_frequency = model$bootstrap_info$cp_frequency
  y_values = cp_frequency*scale_multiplier
  if(pkg == "base"){
    plot(index_values, y_values,
         xlab = index_variable_name,
         ylab = paste0("Detection ", scale),
         main = paste0("Bootstrap estimated detection ",
                       scale, " vs ", index_variable_name))
    lines(index_values, y_values)
    abline(v = model$changepoints, col = "red", lty = "dashed")
  }

}
