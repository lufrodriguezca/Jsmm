# R/number_of_captures.R

#' Calculate the number of captures
#'
#'This function computes captured individuals for each species. This function is used for evaluating the model fit (posterior predictive data).
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @param captures Matrix with the capture data.
#'
#' @return A vector with the computations for each species.
#' @export
#'
number_of_captures = function(m, captures){
  if(is.null(nrow(captures))){
    n.captured = rep(0, m[["ns"]])
  } else {
    n.captured = rep(NA, m[["ns"]])
    captured = m[["releases"]][captures[, 1], 1]
    for(i in 1:m[["ns"]]){
      n.captured[i] = sum(captured == i)
    }
  }
  n.captured
}
