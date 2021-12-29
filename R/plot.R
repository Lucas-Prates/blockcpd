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
    plot(sf, xlim = c(index_values[1], index_values[ncol]), do.points = F, xaxs = "i",
         xlab = index_variable_name, ylab = parameter,
         main = paste("Block plot for", parameter, "parameter"))
  }

}

#' @title
#' Plot to aid choosing regularization constant
#'
#' @description
#' Plots how the number of change points estimated by the given model vary
#' with the regularization constant lambda. The graphic is combined with the
#' "First Long Flat Interval" heuristics to choose a constant value.
#' It is similar to the Elbow method used in clustering. The values of the
#' constant range from 'lambda_left' to 'lambda_right', increasing by 'step'.
#' For each value, the function \link[=fit_blockcpd]{fit_blockcpd} is run
#' with arguments 'blockcpd_args'.
#'
#' @param data_matrix Data frame or matrix containing the data set to be
#' segmented.
#' @param lambda_left Left most value of lambda. Must be non-negative.
#' @param lambda_right Right most value of lambda. Must be non-negative and
#' greater than lambda_left.
#' @param step Value by which lambda will be increased. Must be greater than 0.
#' @param blockcpd_args A list with argument values for the
#' \link[=fit_blockcpd]{fit_blockcpd} function. The list keys must be the
#' arguments names. It must *not* contain the argument 'lambda' or
#' 'data_matrix'.
#' @param pkg Graphical package to be used for plotting. Current values are
#' "base".
#'
#' @return Along with the plot, it returns a list containing the lambda values,
#' number of change points per lambda, the negative log likelihood per lambda
#' and the blockcpd_args.
#' @export
flatplot = function(data_matrix, lambda_left = 0, lambda_right = 10, step = 0.5,
                    blockcpd_args = list(), pkg = "base"){

  # check input
  args_to_check = list(blockcpd_args = blockcpd_args,
                       lambda_left = lambda_left,
                       lambda_right = lambda_right,
                       step = step)
  do.call(check_input, list(caller = as.character(match.call()[[1]]),
                            args_to_check = args_to_check))


  lambda_set = seq(lambda_left, lambda_right, step)
  lambda_set_len = length(lambda_set)
  ncp = numeric(lambda_set_len)
  neg_loglike = numeric(lambda_set_len)
  call_arg = c(list(data_matrix = data_matrix, lambda = lambda_set[1]),
               blockcpd_args)
  call_arg$skip_input_check = TRUE # skips input to avoid error
  for(i in 1:lambda_set_len){
    call_arg$lambda = lambda_set[i]
    model = do.call(fit_blockcpd, call_arg)
    ncp[i] = model$ncp
    neg_loglike[i] = model$neg_loglike
  }
  if(pkg == "base"){
    plot(lambda_set, ncp,
         xlab = "lambda", ylab = "Number of change points",
         main = "Flatplot - number of change points per regularization constant value")
    lines(lambda_set, ncp)
  }
  flatplot_info = list(lambda = lambda_set, ncp = ncp,
                       neg_loglike = neg_loglike)

  invisible(flatplot_info)
}
