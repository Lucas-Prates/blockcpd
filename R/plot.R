#' @title
#' Plot for blockcpd object
#'
#' @description
#' Plots the selected parameters in a blocked fashion.
#'
#' @param x A fitted blockcpd S3 object provided by the
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
#' @param ... Other parameters
#' @rdname plot.blockcpd
#' @export
plot.blockcpd = function(x, ..., parameter = NULL,
                         index_values = NULL, index_variable_name = "Index",
                         pkg = "base"){
  blockcpd_obj = x
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
  do.call(check_input, list(caller = "plot.blockcpd",
                            args_to_check = args_to_check))

  ncp = blockcpd_obj$ncp
  changepoints = blockcpd_obj$changepoints
  parameter_vec = blockcpd_obj$parameters[[parameter]]
  if(pkg == "base"){
    sf = stats::stepfun(index_values[changepoints], parameter_vec)
    stats::plot.stepfun(sf, xlim = c(index_values[1], index_values[ncol]),
                        do.points = F,
                        xaxs = "i", verticals = F,
                        xlab = index_variable_name, ylab = parameter,
                        main = paste("Block plot for", parameter, "parameter"))
    graphics::abline(v = blockcpd_obj$changepoints, col = "red", lty = "dashed")
  }

}

#' @title
#' Plot for graphical selection of the constant
#'
#' @description
#' Plots the output of a frv object. It shows how the number of change-points
#' estimated by the given model vary with the regularization constant lambda.
#' Graphical inspection can be used to choose a proper value for the constant.
#' The suggestion is to pick a value in which the curve starts to "flat-out"
#' @param x An object returned from the function
#' \link[=select_frv]{select_frv}
#' @param pkg Graphical package to be used for plotting. Current values are
#' "base".
#' @param ... Other parameters
#' @rdname plot.frv
#' @export
plot.frv = function(x, ..., pkg = "base"){
  frv_obj = x
  args_to_check = list(frv_obj = frv_obj)
  do.call(check_input, list(caller = "plot.frv",
                            args_to_check = args_to_check))
  lambda_set = frv_obj$lambda_set
  ncp = frv_obj$ncp
  suggested_lambda = frv_obj$suggested_lambda
  suggested_ncp = frv_obj$suggested_ncp
  if(pkg == "base"){
    plot(lambda_set, ncp,
         xlab = "lambda", ylab = "Number of change points",
         main = "Number of change points vs penalization constant")
    graphics::lines(lambda_set, ncp)
    if(!is.null(suggested_lambda)){
      graphics::abline(v = suggested_lambda, col = 'red', lty = "dashed")
      graphics::abline(h = suggested_ncp, col = 'red', lty = "dashed")
    }
  }
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
  do.call(check_input, list(caller = "confidence_plot",
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
    graphics::lines(index_values, y_values)
    graphics::abline(v = model$changepoints, col = "red", lty = "dashed")
  }

}
