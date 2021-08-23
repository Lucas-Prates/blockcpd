# Utility function for rcpd
# Generates NA flags for data sampled in rcpd
#
#
generate_na = function(n_samples, prob){
  sample = stats::rbinom(n_samples, size = 1, prob = 1 - prob)
  sample[sample == 0] = NA
  return(sample)
}

#' @title
#' Sampler for the CPD Block Model
#'
#' @description
#' Creates a \eqn{n \times m} matrix with \eqn{k} change points. In between
#' change points, the random variables are i.i.d. sampled from the given family
#' and parameters
#'
#' @param n Sample size.
#' @param m Array length.
#' @param k Number of change points.  The number of blocks is \eqn{k + 1}. It is
#' overridden if changepoints is non-NULL.
#' @param family The family model to be sampled. The families currently
#' implemented are:
#'
#' \itemize{
#'  \item[bernoulli] Sample independent Bernoullis with probability parameter
#'  of the block
#'  \item[normal] Sample independent Normal with mean and variance specified
#'  by the block.
#'  \item[binaryMarkov] Samples a two state Markov Chain process with transition
#'  matrix defined by the block.
#'  \item[exponential] Sample independent Exponential with scale parameter
#'  defined by the block.
#'  \item[poisson] Sample independent Poisson with rate parameter defined
#'  by the block.
#' }
#' @param parameters List of parameters containing \eqn{k + 1} dimensional
#' parameter vectors of each block. If NULL, the parameters are sampled
#' randomly.
#'
#' @param changepoints A sorted vector of size \eqn{k} containing integers as
#'  change point locations. The change points are between 1 and \eqn{m-1}. If
#'  NULL, the change points are sampled uniformly in \eqn{[1, m-1]}.
#' @param  prob_NA Probability of each entry of being NA. Default is 0.
#'
#' @export
rcpd = function(n = 100,
                m = 50,
                k = 1,
                family = "bernoulli",
                parameters = NULL,
                changepoints = NULL,
                prob_NA = 0) {

  # Setup variables if NULL or check for input errors if the user specified the
  # arguments
  #---------
  if (is.null(changepoints)) {
    if ((k >= m)||(k < 0)){
      stop(paste0("Input error! The number of change points k must be between ",
                  "0 and m."))
    }
    if(family == "binaryMarkov"){
      k = min(k, floor(m/2))
      changepoints = sort(c(0, sample(2:(m-2), k, replace = FALSE), m))
      while(min(diff(changepoints)) == 1){
        changepoints = sort(c(0, sample(2:(m-2), k, replace = FALSE), m))
      }
    }
    else{
      changepoints = sort(c(0, sample(1:(m-1), k, replace = FALSE), m))
    }
  }

  # Check if the change point vector provided is valid
  # and append auxiliary change points
  else {

    if ((any(changepoints <= 0)) || (any(changepoints >= m))) {
      stop("Input error! Change point vector entries must vary from 1 to m-1.")
    }

    k = length(changepoints)
    # Auxiliary change points for sampling
    changepoints = c(0, sort(changepoints), m)

  }


  #---------
  IMPLEMENTED_FAMILIES = c("normal", "bernoulli", "binaryMarkov",
                           "exponential", "poisson")

  if (!(family %in% IMPLEMENTED_FAMILIES)){
    stop(paste0("Input error! The argument 'family' provided is not
                 the list of possible families."))
  }
  if(family == "bernoulli"){
    if (is.null(parameters)) {
      parameters = list(prob = stats::runif(k+1))
    }
    sampler = function(n, i){
      return(stats::rbinom(n , size = 1, parameters[[1]][i]))
    }
  }
  if(family == "normal"){
    if (is.null(parameters)) {
      mean_vec = stats::rnorm(k + 1, 0, 10)
      var_vec = stats::rexp(k + 1)
      parameters = list(mean = mean_vec, var = var_vec)
    }
    sampler = function(n, i){
      return(stats::rnorm(n,
                          parameters[[1]][i],
                          sqrt(parameters[[2]][i]))
      )
    }
  }


  # For this family in specific, the sampler has a different form
  # since it is not independent
  if(family == "binaryMarkov"){
    if (is.null(parameters)) {
      parameters = list(prob00 = stats::runif(k+1),
                        prob11 = stats::runif(k+1))
    }
    data_matrix = rcpd_cpp(family = "binaryMarkov",
                           n,
                           m,
                           changepoints,
                           parameters)
    if(prob_NA != 0){
      data_matrix = data_matrix * generate_na(n*m, prob_NA)
    }
    data_list = list(data_matrix = data.matrix(data_matrix),
                     changepoints = changepoints[-c(1, k+2)], # remove 0 and m
                     parameters = parameters)
    return(data_list)
  }

  if(family == "exponential"){
    if (is.null(parameters)) {
      parameters = list(scale = stats::runif(k+1))
    }
    sampler = function(n, i){
      return(stats::rexp(n, 1.0/parameters[[1]][i]))
    }
  }

  if(family == "poisson"){
    if (is.null(parameters)) {
      parameters = list(rate = stats::runif(k+1))
    }
    sampler = function(n, i){
      return(stats::rpois(n , parameters[[1]][i]))
    }
  }

  # Generate data
  data_matrix = matrix(0, nrow = n, ncol = m)
  for (i in 1:(k + 1)) {
    # A block starts at changepoints[i] + 1 and ends at changepoints[i+1]
    interval = (changepoints[i] + 1):(changepoints[i + 1])
    n_samples = length(interval) * n
    data_matrix[, interval] = sampler(n_samples, i)
    if (prob_NA != 0) {
      data_matrix[, interval] = data_matrix[, interval] * generate_na(n_samples,
                                                                 prob_NA)
    }
  }

  # Put relevant information in a list and remove auxiliary change points
  data_list = list(data_matrix = data.matrix(data_matrix),
                    changepoints = changepoints[-c(1, k+2)], # remove 0 and m
                    parameters = parameters)
  return(data_list)
}
