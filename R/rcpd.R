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
#' Creates a \eqn{nrow \times ncol} matrix with \eqn{ncp} change points.
#' In between change points, the random variables are i.i.d. sampled from the
#' given family and parameters
#'
#' @param nrow Number of rows, or sample size, of the data.
#' @param ncol Number of columns of data matrix. It is the number of variables for
#' each sample.
#' @param ncp Number of change points.  The number of blocks is \eqn{ncp + 1}.
#' It is overridden if changepoints is non-NULL.
#' @param family The family model to be sampled. The families currently
#' implemented are:
#'
#' \itemize{
#'  \item bernoulli: Sample independent Bernoullis with probability parameter
#'  of the block
#'  \item normal: Sample independent Normal with mean and variance specified
#'  by the block.
#'  \item binaryMarkov: Samples a two state Markov Chain process with transition
#'  matrix defined by the block.
#'  \item exponential: Sample independent Exponential with scale parameter
#'  defined by the block.
#'  \item poisson: Sample independent Poisson with rate parameter defined
#'  by the block.
#' }
#' @param parameters List of parameters containing \eqn{ncp + 1} dimensional
#' parameter vectors of each block. If NULL, the parameters are sampled
#' randomly.
#'
#' @param changepoints A sorted vector of size \eqn{ncp} containing integers as
#'  change point locations. The change points are between 1 and \eqn{ncol-1}. If
#'  NULL, the change points are sampled uniformly in \eqn{[1, ncol-1]}.
#' @param  prob_NA Probability of each entry of being NA. Default is 0.
#' @return Returns a list containing 3 elements:
#' #' \itemize{
#'  \item{"data_matrix"} A matrix containing the data.
#'  \item{"changepoints"} A numeric vector containing the change-point locations
#'  \item{"parameters"} A list whose keys are the parameters names and the
#'  values are vectors containing the parameter for each block.
#' }
#' @examples
#' td = rcpd(nrow = 20, ncol = 10) # 20 Bernoulli series of size 10 with 1 change-point
#' td = rcpd(nrow = 10, ncol = 100, ncp = 5,
#'           family = "normal") # 10 normal series of size 100 with 5 change-points
#' td = rcpd(nrow = 1000, ncol = 100, changepoints = c(10, 40, 79)) # choosing change-points locations
#' td = rcpd(nrow = 100, ncol = 15, ncp = 2, family = "normal",
#'           parameters = list(mean = c(1, 2, 3), var = c(4, 5, 6))) # choosing parameters
#' @export
rcpd = function(nrow = 100,
                ncol = 50,
                ncp = 1,
                family = "bernoulli",
                parameters = NULL,
                changepoints = NULL,
                prob_NA = 0) {

  # Setup variables if NULL or check for input errors if the user specified the
  # arguments
  #---------
  IMPLEMENTED_FAMILIES = c("normal", "bernoulli", "binaryMarkov",
                           "exponential", "poisson")

  args_to_check = list(ncp = ncp, ncol = ncol, nrow = nrow,
                       changepoints = changepoints,
                       family = family,
                       IMPLEMENTED_FAMILIES = IMPLEMENTED_FAMILIES)
  do.call(check_input, list(caller = "rcpd", args_to_check = args_to_check))

  if (is.null(changepoints)) {
    if(family == "binaryMarkov"){
      ncp = min(ncp, floor(ncol/2))
      changepoints = sort(c(0, sample(2:(ncol-2), ncp, replace = FALSE), ncol))
      while(min(diff(changepoints)) == 1){
        changepoints = sort(c(0, sample(2:(ncol-2), ncp, replace = FALSE), ncol))
      }
    }
    else{
      changepoints = sort(c(0, sample(1:(ncol-1), ncp, replace = FALSE), ncol))
    }
  }

  # Check if the change point vector provided is valid
  # and append auxiliary change points
  else {
    ncp = length(changepoints)
    # Auxiliary change points for sampling
    changepoints = c(0, sort(changepoints), ncol)
  }


  #---------

  if(family == "bernoulli"){
    if (is.null(parameters)) {
      parameters = list(prob = stats::runif(ncp + 1))
    }
    sampler = function(nrow, i){
      return(stats::rbinom(nrow , size = 1, parameters[[1]][i]))
    }
  }
  if(family == "normal"){
    if (is.null(parameters)) {
      mean_vec = stats::rnorm(ncp + 1, 0, 10)
      var_vec = stats::rexp(ncp + 1)
      parameters = list(mean = mean_vec, var = var_vec)
    }
    sampler = function(nrow, i){
      return(stats::rnorm(nrow,
                          parameters[[1]][i],
                          sqrt(parameters[[2]][i]))
      )
    }
  }


  # For this family in specific, the sampler has a different form
  # since it is not independent
  if(family == "binaryMarkov"){
    if (is.null(parameters)) {
      parameters = list(prob00 = stats::runif(ncp + 1),
                        prob11 = stats::runif(ncp + 1))
    }
    data_matrix = rcpd_cpp(family = "binaryMarkov",
                           nrow,
                           ncol,
                           changepoints,
                           parameters)
    if(prob_NA != 0){
      data_matrix = data_matrix * generate_na(nrow*ncol, prob_NA)
    }
    data_list = list(data_matrix = data.matrix(data_matrix),
                     changepoints = changepoints[-c(1, ncp + 2)], # remove aux cp
                     parameters = parameters)
    return(data_list)
  }

  if(family == "exponential"){
    if (is.null(parameters)) {
      parameters = list(scale = stats::runif(ncp + 1))
    }
    sampler = function(nrow, i){
      return(stats::rexp(nrow, 1.0/parameters[[1]][i]))
    }
  }

  if(family == "poisson"){
    if (is.null(parameters)) {
      parameters = list(rate = stats::runif(ncp + 1))
    }
    sampler = function(nrow, i){
      return(stats::rpois(nrow, parameters[[1]][i]))
    }
  }

  # Generate data
  data_matrix = matrix(0, nrow = nrow, ncol = ncol)
  for (i in 1:(ncp + 1)) {
    # A block starts at changepoints[i] + 1 and ends at changepoints[i+1]
    interval = (changepoints[i] + 1):(changepoints[i + 1])
    n_samples = length(interval) * nrow
    data_matrix[, interval] = sampler(n_samples, i)
    if (prob_NA != 0) {
      data_matrix[, interval] = data_matrix[, interval] * generate_na(n_samples,
                                                                 prob_NA)
    }
  }

  # Put relevant information in a list and remove auxiliary change points
  data_list = list(data_matrix = data.matrix(data_matrix),
                    changepoints = changepoints[-c(1, ncp + 2)], # remove aux cp
                    parameters = parameters)
  return(data_list)
}
