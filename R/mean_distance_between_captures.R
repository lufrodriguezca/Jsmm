# R/mean_distance_between_captures.R

#' Compute the mean distance between captures
#'
#' This function computes the mean euclidean distance between captures for each species. This function is used for evaluating the model fit (posterior predictive data).
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @param captures Matrix with the capture data.
#'
#' @return A vector with the computations for each species.
#' @export
#'
mean_distance_between_captures = function(m, captures){
  if(is.null(nrow(captures))){
    me = rep(NA, m[["ns"]])
  } else {
    ns = m[["ns"]]
    id = captures[, 1]
    sp =  m[["releases"]][id, 1]
    ce = captures[, 2]
    re = m[["releases"]][captures[, 1], 2]
    cl = m[["observation_effort"]][["captures"]][["location"]][ce]
    rl = m[["observation_effort"]][["releases"]][["location"]][re]
    nl = length(m[["observation_effort"]][["locations"]])
    xy = matrix(nrow = nl, ncol = 2)
    for(i in 1:nl){
      sel = m[["site_nodes"]][, i] == 1
      xy[i, ] = colMeans(m[["domain"]][["triangulation"]][["node"]][sel, ])
    }
    dxy = xy[cl, ] - xy[rl, ]
    di = sqrt(dxy[, 1]^2+dxy[, 2]^2)
    me = rep(NA,ns)
    for(i in 1:ns){
      sel = which(sp == i)
      if(length(sel) > 0) me[i] = mean(di[sel])
    }
  }
  me
}
