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
#' @param lambda The penalization constant. Must be a unique non-negative
#' numeric value.
#' @param pen_func Regularization function used for fitting, with default as the
#' BIC. For user specified functions, check the template in the
#' \link[=toy_regularization]{regularization} regularization.rd file.
#' @param min_block_size Minimum block size allowed. Default is 1, and the value
#' must be smaller or equal to ncol.
#' @param max_blocks An integer greater than 0 that specify the maximum number
#' of blocks fitted by the algorithm. It is only used if dynseg is specified in
#' the "method" argument.
#' @param bootstrap A flag to decide if bootstrap computations for the
#' estimation of the probability of each index being detected as a change point.
#' It also provides a sample of all the metrics implemented computed with
#' respect to the final change point set estimated.
#' @param boostrap_samples Number of bootstrap samples.
#' @param skip_input_check Flag indicating if input checking should be skipped.
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
#'  \item{"ncp"} number of change points estimated;
#'  \item{"metadata"} Arguments passed to fit the model;
#'  \item{"bootstrap_info"} if bootstrap argument is true, this contains a list
#'  of the metrics for each bootstrap sample, and contains the estimated
#'  probability of each index being detected as a change point;
#' }
#' @export
fit_blockcpd = function(data_matrix,
                        method = "hierseg",
                        family = "bernoulli",
                        lambda = 1.0,
                        pen_func = bic_loss,
                        min_block_size = 1L,
                        max_blocks = NULL,
                        bootstrap = FALSE,
                        bootstrap_samples = 100L,
                        bootstrap_progress = FALSE,
                        skip_input_check = FALSE) {

  nrow = nrow(data_matrix)
  ncol = ncol(data_matrix)
  # input check
  if(!skip_input_check){
    args_to_check = list(method = method,
                         family = family,
                         ncol = ncol,
                         min_block_size = min_block_size,
                         lambda = lambda)
    do.call(check_input, list(caller = "fit_blockcpd",
                              args_to_check = args_to_check))
  }
  if(is.null(max_blocks)){max_blocks = ncol}

  methodcall_name = paste0("compute_", method) # method name in package
  suff_stats = compute_suff_stats(data_matrix, family)

  # Fits model for the given lambda
  fit_arguments   = list(suff_stats = suff_stats,
                         family = family,
                         nrow = nrow,
                         ncol = ncol,
                         lambda = lambda,
                         pen_func = pen_func,
                         min_block_size = min_block_size,
                         max_blocks = max_blocks)

  model = do.call(methodcall_name, fit_arguments)

  ### ------------------------------------------------- ###
  # Bootstrap computation
  # ??? Separate this in a new function ???
  if(bootstrap){
    cp_frequency = rep(0, ncol) # Frequency of an index being detected as a change point
    haus_boot = rep(0, bootstrap_samples)
    symdiff_values = rep(0, bootstrap_samples)
    randindex_values = rep(0, bootstrap_samples)
    jaccard_values = rep(0, bootstrap_samples)
    ncp_values = rep(0, bootstrap_samples)

    for(i in 1:bootstrap_samples){
      boot_samp = sample(nrow, nrow, replace = TRUE)
      suff_stats_boot = compute_suff_stats(data_matrix[boot_samp, ], family)
      fit_arguments_boot = list(suff_stats = suff_stats_boot,
                                family = family,
                                nrow = nrow,
                                ncol = ncol,
                                lambda = lambda,
                                pen_func = pen_func,
                                min_block_size = min_block_size,
                                max_blocks = max_blocks)
      model_boot = do.call(methodcall_name, fit_arguments_boot)
      # For each index (column) detected as change point, increment it
      cp_frequency[model_boot$changepoints] = cp_frequency[model_boot$changepoints] + 1

      # Metrics evaluation with respect to final estimate of change point set
      haus_boot[i] = compute_hausdorff(model$changepoints,
                                       model_boot$changepoints)
      symdiff_values[i] = compute_symdiff(model$changepoints,
                                        model_boot$changepoints)
      randindex_values[i] = compute_rand(model$changepoints,
                                  model_boot$changepoints, ncol)
      jaccard_values[i] = compute_jaccard(model$changepoints,
                                        model_boot$changepoints)
      ncp_values[i] = length(model_boot$changepoints)
      if(bootstrap_progress){
        cat(paste0("\r Iteration ", i, " of ", bootstrap_samples))
      }
    }
    if(bootstrap_progress){cat("\n")}

    bootstrap_info = list(bootstrap_samples = bootstrap_samples,
                          cp_frequency = cp_frequency,
                          symdiff_values = symdiff_values,
                          randindex_values = randindex_values,
                          jaccard_values = jaccard_values,
                          ncp_values = ncp_values)

    model$bootstrap_info = bootstrap_info
  }else{model$bootstrap_info = NULL}
  ### ------------------------------------------------- ###

  # Wrap results in a S3 class called blockcpd
  model$metadata = list(method = method,
                        family = family,
                        nrow = nrow,
                        ncol = ncol,
                        bootstrap = bootstrap,
                        lambda = lambda,
                        pen_func = pen_func,
                        min_block_size = min_block_size,
                        max_blocks = max_blocks)
  class(model) = "blockcpd"

  return(model)

}
