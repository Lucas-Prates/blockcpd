#' @title
#' Estimates crudely the curvature of the graph of a function
#'
#' @description Utility method used in elbow_plot. Given an interval
#' and function values, estimates the curvature of a function at the
#' interior points. The domain values are apart by a distance of 'step'.
#' @noRd
compute_max_curvature = function(dom_values, fvalues, step){
  max_curvature = -Inf
  len_domain = length(dom_values)
  for(i in 2:(len_domain-1)){
    if((2 < i)&&(i < len_domain-1)){
      ncp_diff1 = (fvalues[i+1] - fvalues[i-1])/(2*step)# estimating first derivative at i
      ncp_diff2 = (fvalues[i+2] - fvalues[i-2])/(4*step*step)# estimating second derivative at i
    }
    if(i == 2){
      ncp_diff1 = (fvalues[i+1] - fvalues[i-1])/(2*step)# estimating first derivative at i
      ncp_diff2 = (fvalues[i+2] - fvalues[i-1])/(2*step*step)# estimating second derivative at i
    }
    else{
      ncp_diff1 = (fvalues[i+1] - fvalues[i-1])/(2*step)# estimating first derivative at i
      ncp_diff2 = (fvalues[i+1] - fvalues[i-2])/(2*step*step)# estimating second derivative at i
    }
    curvature = abs(ncp_diff2)/((1 + ncp_diff1*ncp_diff1)^(3/2))
    if(max_curvature < curvature){
      max_curvature = curvature
      argmax_curv = dom_values[i]
      fmax_curv = fvalues[i]
    }
  }
  curvature_info = list(argmax_curv = argmax_curv,
                        fmax_curv = fmax_curv)
  return(curvature_info)
}
