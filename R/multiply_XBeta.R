# R/multiply_XBeta.R

#' @noRd
#'
multiply_XBeta = function(m, Beta, i){
  X = m[["X"]]
  included = (m$nc > 0)
  vars = names(included)
  nele =  nrow(m[["domain"]][["triangulation"]][["ele"]])

  XBeta = list()
  for(j in 1:5){
    if(included[vars[[j]]]){
      if(m[["dynamic_covs_movement"]]){
        XBeta[[vars[[j]]]] = X[[vars[[j]]]][[i]] %*% t(matrix(Beta[[vars[[j]]]], nrow = m[["ns"]]))
      } else {
        XBeta[[vars[[j]]]] = X[[vars[[j]]]] %*% t(matrix(Beta[[vars[[j]]]], nrow = m[["ns"]]))
      }
    }
  }

  if(included[["diffusion"]]) XBeta[["diffusion"]] = exp(XBeta[["diffusion"]])
  if(included[["mortality"]]) XBeta[["mortality"]] = exp(XBeta[["mortality"]])
  if(included[["habitat_preference"]]) XBeta[["habitat_preference"]] = exp(XBeta[["habitat_preference"]])
  if(is.null(XBeta[["habitat_preference"]])) XBeta[["habitat_preference"]] = matrix(1, nrow = 3*nele, ncol = m[["ns"]])

  if(included[["observation"]] && m[["observation_effort"]][["method"]] == "CCP"){
    tmp = matrix(0, ncol = m[["ns"]], nrow = 3*nele)
    ne = length(m[["observation_effort"]][["captures"]][["location"]]) #number of observation events
    for(j in seq_len(ne)){
      if(m[["dynamic_covs_observation"]]){
        o_times = X[["observation_times"]][[j]]
      } else {
        o_times =  X[["observation_times"]][j, ]
      }
      if(i >= o_times[1] & i< o_times[2]){
        if(m[["dynamic_covs_observation"]]){
          x = X[["observation"]][[j]]
        } else {
          x = X[["observation"]][j, ]
        }

        r = x%*%t(Beta[["observation"]])

        site = m[["observation_effort"]][["captures"]][["location"]][j]
        elements = m[["observation_effort"]][["locations"]][[site]]
        z.one = c(elements,elements + nele,elements + 2*nele)
        for(sp in 1:m[["ns"]]) tmp[z.one, sp] = tmp[z.one, sp] + exp(r[sp])
      }
    }

    XBeta[["observation"]] = tmp
  }

  return(XBeta)
}

