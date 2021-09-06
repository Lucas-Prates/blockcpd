#' @title
#' Fits a blockcpd model
#'
#' @description
#' Fits a blockcpd model to find the best segmentation of the data into blocks.
#' Variables in each block have the same distribution and parameter,
#' and consecutive blocks have different parameters.
#'
#' @param data_matrix Data frame or matrix containing the data set to be
#' segmented. There is no verification if the entries correspond to the model
#' specified by the "family" argument, such as entries different than 0, 1 or NA
#' for the bernoulli family.
#' @param method The method that will be used to fit the model. The current
#' implemented models are:
#' \itemize{
#'  \item[\link[=compute_hierseg]{hierseg}] Hierarchical segmentation, also
#'  known as binary segmentation;
#'  \item[\link[=compute_dynseg]{dynseg}] Dynamical programming segmentation.
#' }
#' @param family The name of the family to detect changes in parameters. Should
#' be passed as a string. The families  currently implemented are:
#'
#' \itemize{
#'  \item{"bernoulli"} The model assumes that data comes from a Bernoulli
#'  distribution. For each block, the algorithm estimates the probability
#'  paramater. Each entry should be binary.
#'  \item{"normal"} The model assumes data comes fro ma Normal distribution with
#'  unknown mean and variance. For each block, the algorithms estimates the
#'  mean and variance parameter. Each entry should be numeric.
#'  \item{"binaryMarkov"} The model assumes that data comes from two states (0, 1)
#'  Markov Chain. For each block, the algorithm estimates the 2x2 transition
#'  matrix. Each entry should be binary. At the boundary of the blocks, the
#'  transition is defined using the parameters of the next (new) block. For
#'  instance, consider a block defined from a to c, followed a block from c + 1
#'  to b (including the extremes). By definition, c is a change point, and the
#'  transition from X_c to X_{c} + 1 is defined by the parameters on c + 1 to
#'  b.
#'  \item{"exponential"} The model assumes that data comes from an Exponential
#'  distribution. For each block, the algorithm estimates the scale parameter,
#'  that is, the inverse of the rate. Each entry should be numeric and positive.
#'  \item{"poisson"} The model assumes that data comes from a Poisson distribution
#'  For each block, the algorithm estimates the rate paramater. Each entry
#'  should an positive integer.
#' }
#'
#' @param lambda The penalization constant. A list of penalization constant can
#' be passed if the argument "select_lambda" is TRUE.
#' @param select_lambda A flag to decide if the BIC criterion must be used to
#' choose the best lambda from the list of values provided in the "lambda"
#' argument.
#' @param pen_func Regularization function used for fitting, with default as the
#' BIC. For user specified functions, check the template in the
#' \link[=toy_regularization]{regularization} regularization.rd file.
#' @param max_blocks An integer greater than 0 that specify the maximum number
#' of blocks fitted by the algorithm. It is only used if dynseg is specified in
#' the "method" argument.
#' @param bootstrap A flag to decide if bootstrap computations for the
#' estimation of the probability of each index being detected as a change point.
#' It also provides a sample of all the metrics implemented computed with
#' respect to the final change point set estimated.
#' @param boostrap_rep Number of bootstrap repetitions.
#'
#' @return The function returns a S3 object of the type blockcpd.
#' \itemize{
#'  \item{"changepoints"} a list containing the set of estimated change points;
#'  \item{"parameters"} a list containing the estimated parameters for each
#'  block. In the case of multiple parameters, it provides a list of lists,
#'  where each sub list refers to the parameter that names the list;
#'  \item{"loss"} the final loss evaluated on the entire data set for the
#'  returned model;
#'  \item{"neg_loglike"} The negative log likelihood of the model;
#'  \item{"n_cp"} number of change points estimated;
#'  \item{"metadata"} Arguments passed to fit the model;
#'  \item{"bootstrap_info"} if bootstrap argument is true, this contains a list
#'  of the metrics for each bootstrap sample, and contains the estimated
#'  probability of each index being detected as a change point;
#'  \item{"lambda_info"} If select_lambda argument is true, this contains a list
#'  summarizing the the BIC values and the best value for lambda.
#' }
#' @export
fit_blockcpd = function(data_matrix,
                        method = "hierseg",
                        family = "bernoulli",
                        lambda = 1,
                        select_lambda = FALSE,
                        pen_func = bic_loss,
                        max_blocks = NULL,
                        bootstrap = FALSE,
                        bootstrap_rep = 100L,
                        bootstrap_progress = FALSE) {

  ### Check inputs
  IMPLEMENTED_METHODS = c("hierseg", "dynseg")

  IMPLEMENTED_FAMILIES = c("normal", "bernoulli", "binaryMarkov",
                           "exponential", "poisson")

  if ( !(method %in% IMPLEMENTED_METHODS) ) {
    stop("Error! The 'method' argument provided is not implemented!")
  }

  if ( !(family %in% IMPLEMENTED_FAMILIES) ) {
    stop("Error! The 'family' argument provided is not implemented!")
  }
  ###

  n = nrow(data_matrix)
  m = ncol(data_matrix)
  if(is.null(max_blocks)){max_blocks = m}

  methodcall_name = paste0("compute_", method) # method name in package
  suff_stats = compute_suff_stats(data_matrix, family)

  # Fits a model for each lambda and chooses the best using the BIC criterion
  if(select_lambda){
    bic_value = rep(0, length(lambda))
    best_bic = Inf
    best_i = lambda[1]
    for(i in 1:length(lambda)){
      lambda_i = lambda[i]
      fit_arguments = list(suff_stats = suff_stats,
                           family = family,
                           n = n,
                           m = m,
                           lambda = lambda_i,
                           pen_func = pen_func,
                           max_blocks = max_blocks)

      model_lambda = do.call(methodcall_name, fit_arguments)
      n_param = model_lambda$n_cp + 1
      bic_value[i] = 2*model_lambda$neg_loglike + log(n)*n_param
      # Avoid recomputation of the model
      if(bic_value[i] < best_bic){
        model = model_lambda
        best_lambda = lambda[i]
        best_bic = bic_value[i]
      }
    }
  }
  else{
    # Fits model for the unique given lambda
    # If a list of lambdas is passed, only the first value is considered
    best_lambda = lambda # for bootstrap argument consistency
    fit_arguments   = list(suff_stats = suff_stats,
                           family = family,
                           n = n,
                           m = m,
                           lambda = lambda,
                           pen_func = pen_func,
                           max_blocks = max_blocks)

    model = do.call(methodcall_name, fit_arguments)
  }


  ### ------------------------------------------------- ###
  # Bootstrap computation
  # Notice that best_lambda is used in the bootstrap computations
  # ??? Separate this in a new function ???
  if(bootstrap){
    cp_freq = rep(0, m) # Frequency of an index being detected as a change point
    haus_boot = rep(0, bootstrap_rep)
    symdiff_boot = rep(0, bootstrap_rep)
    rand_boot = rep(0, bootstrap_rep)
    jaccard_boot = rep(0, bootstrap_rep)
    ncp_boot = rep(0, bootstrap_rep)

    for(i in 1:bootstrap_rep){
      boot_samp = sample(n, n, replace = TRUE)
      suff_stats_boot = compute_suff_stats(data_matrix[boot_samp, ], family)
      fit_arguments_boot = list(suff_stats = suff_stats_boot,
                                family = family,
                                n = n,
                                m = m,
                                lambda = best_lambda,
                                pen_func = pen_func,
                                max_blocks = max_blocks)
      model_boot = do.call(methodcall_name, fit_arguments_boot)
      # For each index (column) detected as change point, increment it
      cp_freq[model_boot$changepoints] = cp_freq[model_boot$changepoints] + 1

      # Metrics evaluation with respect to final estimate of change point set
      haus_boot[i] = compute_hausdorff(model$changepoints,
                                       model_boot$changepoints)
      symdiff_boot[i] = compute_symdiff(model$changepoints,
                                        model_boot$changepoints)
      rand_boot[i] = compute_rand(model$changepoints,
                                  model_boot$changepoints, m)
      jaccard_boot[i] = compute_jaccard(model$changepoints,
                                        model_boot$changepoints)
      ncp_boot[i] = length(model_boot$changepoints)
      if(bootstrap_progress){
        cat(paste0("\r Iteration ", i, " of ", bootstrap_rep))
      }
    }
    if(bootstrap_progress){cat("\n")}
    bootstrap_info = list(b_samples = bootstrap_rep,
                           cp_freq = cp_freq,
                           symdiff_boot = symdiff_boot,
                           rand_boot = rand_boot,
                           jaccard_boot = jaccard_boot,
                           ncp_boot = ncp_boot)

    model$bootstrap_info = bootstrap_info
  }else{model$bootstrap_info = NULL}
  ### ------------------------------------------------- ###

  # Wrap results in a S3 class called blockcpd
  model$metadata = list(method = method,
                        family = family,
                        n = n,
                        m = m,
                        bootstrap = bootstrap,
                        lambda = lambda,
                        pen_func = pen_func,
                        max_blocks = max_blocks)
  if(select_lambda){
    model$lambda_info = list(lambda = lambda,
                             bic_value = bic_value,
                             best_lambda = best_lambda)
  }
  class(model) = "blockcpd"

  return(model)

}
