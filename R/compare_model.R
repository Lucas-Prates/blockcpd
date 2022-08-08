#' @title
#' Compare or evaluate model performance with respect
#' to other model or ground truth
#'
#' @description Compares or evaluates model estimated change point set against
#' another model or ground truth. The comparison is made using common metrics
#' to compare clusters. The metrics provided are
#'\itemize{
#'  \item "hausdorff": Hausdorff Distance metric;
#'  \item "rand": Rand Index ;
#'  \item "symdiff": Symmetric difference metric;
#'  \item "jaccard": Jaccard similarity index.
#' }
#'
#' @param model1 The first blockcpd object or list of sorted integers
#' representing the change point set.
#' @param model2 The second blockcpd object or list of sorted integers
#' representing the change point set.
#' @param ncol The number of variables which the model was fitted on. Only needs
#' to be passed if both arguments are change point sets instead of a blockcpd
#' object.
#' @return Returns a list containing four metrics:
#' \itemize{
#'  \item{"haus"} Hausdorff distance;
#'  \item{"rand"} Rand index;
#'  \item{"symdiff"} Symmetric difference;
#'  \item{"jaccard"} Jaccard index.
#' }
#' @examples
#' model1 = fit_blockcpd(c(0, 1, 0, 1), lambda = 0)
#' model2 = fit_blockcpd(c(0, 1, 0, 1), lambda = Inf)
#' comparison = compare_model(model1, model2)
#' # change-point sets can also be passed directly with ncol
#' compare_model(c(1,2,3,4), c(3), 10)
#' @export

compare_model = function(model1,
                         model2,
                         ncol = NULL){

  blockcpd_flag = FALSE # flag for blockcpd object

  # Model 1 setup
  if (inherits(model1 ,"blockcpd")){
    ncol = model1$metadata$ncol
    cp1 = model1$changepoints
    blockcpd_flag = TRUE
  } else {cp1 = model1} # cp1 must be a change point set

  # Model 2 setup
  if (inherits(model2, "blockcpd")){
    ncol = model2$metadata$ncol
    cp2 = model2$changepoints
    blockcpd_flag = TRUE
  } else {cp2 = model2} # cp2 must be a change point set

  args_to_check = list(ncol = ncol,
                       blockcpd_flag = blockcpd_flag)
  check_input("compare_model", args_to_check)
  # Metrics computations in c++
  cp1 = c(cp1, ncol)
  cp2 = c(cp2, ncol)
  rand = compute_rand(cp1, cp2, ncol)
  haus = compute_hausdorff(cp1, cp2)
  symdiff = compute_symdiff(cp1, cp2)
  jaccard = compute_jaccard(cp1, cp2)

  metrics_list = list(haus = haus,
                      rand = rand,
                      symdiff = symdiff,
                      jaccard = jaccard)

  return(metrics_list)

}
