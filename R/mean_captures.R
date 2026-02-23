# R/mean_captures.R

#' Compute the captured individuals mean
#'
#' This function computes the mean of the captured individual for each species. This function is used for evaluating the model fit (posterior predictive data).
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @param captures Matrix with the capture data.
#'
#' @return A vector with the estimates for each species.
#' @export
#'
mean_captures = function(m, captures){
  n.released = rep(NA, m[["ns"]])
  released = m[["releases"]][, 1]

  for(i in 1:m[["ns"]]){
    n.released[i] = sum(released == i)
  }
  if(is.null(nrow(captures))){
    n.captured = rep(0,m[["ns"]])
  } else {
    n.captured = rep(NA,m[["ns"]])
    captured = m[["releases"]][captures[, 1], 1]
    for(i in 1:m[["ns"]]){
      n.captured[i] = sum(captured == i)
    }
  }
  n.captured/n.released
}


