#' @title
#' Plot for blockcpd object
#'
#' @description
#' Plots the selected parameters in a blocked fashion.
#'
#' @param blockcpd_obj A fitted blockcpd S3 object provided by the
#' \item[\link[=fit_blockcpd]{fit_blockcpd}] function.
#' @param parameter The parameter of the family for which to plot the blocked
#' @export
plot.blockcpd = function(blockcpd_obj, parameter = NULL,
                         library = "base"){
  # check if parameter argument is in family parameter list
  if(is.null(parameter)){
    parameter = names(blockcpd_obj$parameters)[1]
  }
  if(!(parameter %in% names(blockcpd_obj$parameters))){
    stop("Input error! The 'parameter' argument is not a parameter of the family fitted for the blockcpd object!")
  }
  ncol = blockcpd_obj$metadata$ncol
  ncp = blockcpd_obj$ncp
  changepoints = blockcpd_obj$changepoints
  parameter_vec = blockcpd_obj$parameters[[parameter]]
  if(library == "base"){
    sf = stepfun(changepoints, parameter_vec)
    plot(sf, xlim = c(1, ncol), do.points = F, xaxs = "i",
         xlab = "Index", ylab = parameter,
         main = paste("Block plot for", parameter, "parameter"))
  }

}
