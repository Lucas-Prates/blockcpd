#' @title
#' Provides sufficient statistics to fit the model for the given family
#'
#' @description Provides a list of the sufficient statistics for implemented
#' models. For most of them, the cumulative sum of the sum of columns is
#' sufficient.
#'
#' @param data_matrix The data set for which we compute sufficient statistics
#' @param family The name of the family used to fit the model
#' @noRd
compute_suff_stats = function(data_matrix, family){

  if(family == "bernoulli"){
    suff_stats = list(cumsum(colSums(data_matrix, na.rm = T)),
                      cumsum(colSums(!is.na(data_matrix)))
                     )
  }

  if(family == "normal"){
    suff_stats = list(cumsum(colSums(data_matrix, na.rm = T)),
                      cumsum(colSums(data_matrix^2, na.rm = T)),
                      cumsum(colSums(!is.na(data_matrix)))
                     )
  }

  return(suff_stats)
}
