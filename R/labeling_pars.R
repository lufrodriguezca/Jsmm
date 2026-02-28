# R/labeling_pars.R

#' @noRd
#'
labeling_pars = function(m){

  included = rep(NA, 6)
  names(included) = c("diffusion", "advection_1", "advection_2", "mortality", "habitat_preference", "observation")

  model_pars = names(m[["X"]])
  del = c(grep("times", model_pars), grep("trait", model_pars), grep("observation", model_pars))

  if(length(del) > 0){
    model_pars = model_pars[-del]
  }

  parNames = c()
  for(n in model_pars){
    if(m[["dynamic_covs_movement"]]){
      parNames = c(parNames, paste0(n, "_", colnames(m[["X"]][[n]][[1]])))
    }else{
      parNames = c(parNames, paste0(n, "_", colnames(m[["X"]][[n]])))
    }
  }

  if("observation"%in%names(m$X)){
    if(m[["dynamic_covs_observation"]]){
      parNames = c(parNames, paste0(n, "_", colnames(m[["X"]][["observation"]][[1]])))
    }else{
      parNames = c(parNames, paste0("observation_", colnames(m[["X"]][["observation"]])))
    }
  }
  return(parNames)
}
