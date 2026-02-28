# R/keep_used_variables.R

#' @noRd
#'
keep_used_variables = function(cov_data, tr_data, model_formula){

  included.movement.vars = unique(c(all.vars(model_formula[["diffusion"]]),all.vars(model_formula[["advection_1"]]),
                                    all.vars(model_formula[["advection_2"]]),
                                    all.vars(model_formula[["mortality"]]),all.vars(model_formula[["habitat_preference"]])))
  included.observation.vars = all.vars(model_formula[["observation"]])
  included.trait.vars = all.vars(model_formula[["traits"]])

  tr_data = keep_selected_variables(da = tr_data, vars = included.trait.vars)

  cov_data[["n_data"]] = keep_selected_variables(da = cov_data[["n_data"]], vars = included.movement.vars)

  cov_data[["e_data"]] = keep_selected_variables(da = cov_data[["e_data"]], vars = included.movement.vars)

  cov_data[["t_data"]] = keep_selected_variables(da = cov_data[["t_data"]], vars = included.movement.vars)

  if(is.null(cov_data[["t_data"]])) cov_data[["t_data_times"]] = NULL

  if(!is.null(cov_data[["e_t_data"]])){
    tmp = keep_selected_variables(da = cov_data[["e_t_data"]][[1]], vars = included.movement.vars)
    if(is.null(tmp)){
      cov_data[["e_t_data"]] = NULL
      cov_data[["e_t_data_times"]] = NULL
    } else {
      for(i in seq_len(length(cov_data[["e_t_data"]]))){
        cov_data[["e_t_data"]][[i]] = keep_selected_variables(da = cov_data[["e_t_data"]][[i]], vars = included.movement.vars)
      }
    }
  }

  if(!is.null(cov_data[["n_t_data"]])){
    tmp = keep_selected_variables(da = cov_data[["n_t_data"]][[1]], vars = included.movement.vars)
    if(is.null(tmp)){
      cov_data[["n_t_data"]] = NULL
      cov_data[["n_t_data_times"]] = NULL
    }else{
      for(i in seq_len(length(cov_data[["n_t_data"]]))){
        cov_data[["n_t_data"]][[i]] = keep_selected_variables(da = cov_data[["n_t_data"]][[i]], vars = included.movement.vars)
      }
    }
  }

  cov_data[["obs_data"]] = keep_selected_variables(da = cov_data[["obs_data"]], vars = included.observation.vars)

  if(!is.null(cov_data[["obs_t_data"]])){
    for(i in seq_len(length(cov_data[["obs_t_data"]]))){
      cov_data[["obs_t_data"]][[i]][["covs"]] = keep_selected_variables(da = cov_data[["obs_t_data"]][[i]][["covs"]], vars = included.observation.vars)
    }
    if(is.null(cov_data[["obs_t_data"]][[1]][["covs"]])) cov_data[["obs_t_data"]] = NULL
  }

  res = list()
  res[["cov_data"]] = cov_data
  res[["tr_data"]]  = tr_data

  return(res)
}
