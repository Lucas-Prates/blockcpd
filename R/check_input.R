#' @title
#' Checks input from fit_blockcpd
#'
#' @description
#' Checks if the input from \item[\link[=fit_blockcpd]{fit_blockcpd}] passes
#' a sanity check.
check_input = function(method, family, lambda){

  IMPLEMENTED_METHODS = c("hierseg", "dynseg")

  IMPLEMENTED_FAMILIES = c("normal", "bernoulli", "binaryMarkov",
                           "exponential", "poisson")

  if ( !(method %in% IMPLEMENTED_METHODS) ) {
    stop("Input error! The 'method' argument provided is not implemented!")
  }

  if ( !(family %in% IMPLEMENTED_FAMILIES) ) {
    stop("Input error! The 'family' argument provided is not implemented!")
  }

  if ((!is.numeric(lambda))||(length(lambda) != 1)){
    stop("Input error! The 'lambda' argument must be a unique numeric value!")
  }
  if (lambda < 0){
    stop("Input error! The 'lambda' argument must be non-negative!")
  }
}
