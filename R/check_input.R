#' @title
#' Checks input from caller
#'
#' @description
#' Performs a sanity check on the inputs from caller. It stops execution and
#' outputs an error message if arguments are not in conformity with caller
#' method.
check_input = function(caller, args_to_check){

  # checks fit_blockcpd
  #---
  if(caller == "fit_blockcpd"){
    method = args_to_check$method
    family = args_to_check$family
    ncol = args_to_check$ncol
    min_block_size = args_to_check$min_block_size
    lambda = args_to_check$lambda
    IMPLEMENTED_METHODS = c("hierseg", "dynseg")

    IMPLEMENTED_FAMILIES = c("normal", "bernoulli", "binaryMarkov",
                             "exponential", "poisson")

    if(!(method %in% IMPLEMENTED_METHODS) ) {
      stop("Input error! The 'method' argument provided is not implemented!")
    }

    if(!(family %in% IMPLEMENTED_FAMILIES) ) {
      stop("Input error! The 'family' argument provided is not implemented!")
    }

    if((!is.numeric(lambda))||(length(lambda) != 1)){
      stop("Input error! The 'lambda' argument must be a unique numeric value!")
    }
    if((min_block_size > ncol)||(min_block_size <= 0)){
      stop("Input error! The 'min_block_size' argumenst ranges from 1 to ncol!")
    }
    if (lambda < 0){
      stop("Input error! The 'lambda' argument must be non-negative!")
    }
  }
  #---


  # checks plot.blockcpd
  #---
  if(caller == "plot.blockcpd"){
    parameter = args_to_check$parameter
    family_parameters = args_to_check$family_parameters
    is_index_values_numeric = args_to_check$is_index_values_numeric
    length_index_values = args_to_check$length_index_values
    ncol = args_to_check$ncol
    if(!(parameter %in% family_parameters)){
      stop("Input error! The 'parameter' argument is not a parameter of the family fitted for the blockcpd object!")
    }
    if(!is_index_values_numeric){
      stop("Input error! The 'index_values' argument is not a numeric vector!")
    }
    if(length_index_values != ncol){
      stop("Input error! The 'index_values' argument size differs from ncol from the fitted model!")
    }
  }
  #---

  # checks select_frv
  #---
  if(caller == "select_frv"){
    model_args = args_to_check$model_args
    lambda_left = args_to_check$lambda_left
    lambda_right = args_to_check$lambda_right
    step = args_to_check$step
    # checks if model_args is a list
    if(!is.list(model_args)){
      stop("Input error! The 'model_args' argument must be a list!")
    }
    # check if lambda is not in argument list
    if(("lambda" %in% names(model_args))||("data_matrix" %in% names(model_args))){
      stop("Input error! The 'model_args' argument must not contain the 'lambda' or 'data_matrix' as a key!")
    }
    # sanity check on lambda_left, lambda_right
    if((lambda_left >= lambda_right)||(lambda_left < 0)){
      stop("Input error! We must have 0 < 'lambda_left' < 'lambda_right'!")
    }
    # sanity check on step
    if(step <= 0){
      stop("Input error! We must have 'step' > 0")
    }
  }
  #---
  # checks frv_plot
  if(caller == "plot.frv"){
    frv_obj = args_to_check$frv_obj
    if(class(frv_obj) != "frv"){
      stop("Input error! The argument 'frv_obj' must be a frv object!")
    }
  }

  # checks confidence_plot
  if(caller == "confidence_plot"){
    model = args_to_check$model
    scale = args_to_check$scale
    is_index_values_numeric = args_to_check$is_index_values_numeric
    length_index_values = args_to_check$length_index_values
    ncol = args_to_check$ncol
    if(class(model) != "blockcpd"){
      stop("Input error! The argument 'model' must be blockcpd object!")
    }
    if(!model$metadata$bootstrap){
      stop("Input error! Fit the model using 'bootstrap = TRUE'!")
    }
    if(!(scale %in% c("percentage", "probability", "frequency"))){
      stop("Input error! The argument 'scale' must be one of 'percentage', 'probability' or 'frequency'!")
    }
    if(!is_index_values_numeric){
      stop("Input error! The 'index_values' argument is not a numeric vector!")
    }
    if(length_index_values != ncol){
      stop("Input error! The 'index_values' argument size differs from ncol from the fitted model!")
    }
  }
  #---

  # checks rcpd
  if(caller == "rcpd"){
    ncp = args_to_check$ncp
    ncol = args_to_check$ncol
    nrow = args_to_check$nrow
    family = args_to_check$family
    changepoints = args_to_check$changepoints
    IMPLEMENTED_FAMILIES = args_to_check$IMPLEMENTED_FAMILIES

    if (!(family %in% IMPLEMENTED_FAMILIES)){
      stop(paste0("Input error! The argument 'family' provided is not
                   the list of possible families."))
    }
    if ((ncp >= ncol)||(ncp < 0)){
      stop(paste0("Input error! The number of change points ncp must be between ",
                  "0 and ncol-1."))
    }
    if ((any(changepoints <= 0)) || (any(changepoints >= ncol))) {
      stop("Input error! Change point vector entries must vary from 1 to ncol-1.")
    }
  }
  #---
}
