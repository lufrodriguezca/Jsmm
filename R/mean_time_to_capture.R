# R/mean_time_to_capture.R

#'Compute the mean time between captures
#'
#' This function computes the mean time between captures for each species. This function is used for evaluating the model fit (posterior predictive data).
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @param captures Matrix with the capture data.
#'
#' @return A vector with the computations for each species.
#' @export
#'
mean_time_to_capture = function(m,captures){
  if(is.null(nrow(captures))){
    me = rep(NA, m[["ns"]])
  } else {
    ns = m[["ns"]]
    id = captures[, 1]
    sp =  m[["releases"]][id, 1]
    ce = captures[, 2]
    re = m[["releases"]][captures[, 1], 2]
    if(m[["observation_effort"]][["method"]] == "CCP"){
      ct = m[["times"]][m[["observation_effort"]][["captures"]][["time"]][ce, 2]]
    } else{
      ct = m[["times"]][m[["observation_effort"]][["captures"]][["time"]][ce]]
    }
    rt = m[["times"]][m[["observation_effort"]][["releases"]][["time"]][re]]
    me = rep(NA, ns)
    for(i in 1:ns){
      sel = which(sp == i)
      if(length(sel) > 0) me[i] = mean((ct - rt)[sel])
    }
  }
  me
}
