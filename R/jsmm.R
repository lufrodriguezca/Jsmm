# R/jsmm.R

#' Define the Jsmm model
#'
#' Creates a named list with all the necessary components for running the Jsmm analyses.
#'
#' @param domain named list with the elements domain and triangulation. The former is a polygonal object from the class sf, and the latter provides information about the domain triangulation.
#' @param observation_effort named list with spatial and temporal information about releases and captures in the experiment.
#' @param releases named list with information about where and when releases occurred.
#' @param captures named list with information about where and when captures occurred.
#' @param cov_data named list with (some combination of) spatial, temporal or/and spatio-temporal data.
#' @param ns integer representing the number of species included in the model.
#' @param sp_names string vector parameter assigning the species names in the model (optional).
#' @param C phylogenetic relationships matrix between species.
#' @param tr_data dataframe with trait data.
#' @param model_formula named list that consists of (some combination of) the objects diffusion, advection, mortality, habitat_preference, observation, and traits. These follow the R-notation for the model formulae.
#'
#' @return A named list with all the information for running Jsmm analyses.
#' @export
#'
jsmm = function(domain, observation_effort, releases = NULL, captures = NULL,
                cov_data, ns, sp_names = NULL, C = NULL, tr_data = NULL, model_formula){

  if(observation_effort[["method"]] == "ICP"){
    cov_data[["obs_t_data"]] = NULL
  }
  res = keep_used_variables(cov_data = cov_data, tr_data = tr_data, model_formula = model_formula)
  cov_data = res[["cov_data"]]
  tr_data = res[["tr_data"]]

  times = sort(unique(unlist(c(cov_data[["t_data_times"]], cov_data[["n_t_data_times"]], cov_data[["e_t_data_times"]],
                      unlist(observation_effort[["captures"]][["time"]]), observation_effort[["releases"]][["time"]]))))
  if(!is.null(cov_data[["obs_t_data"]])){
    times2 = c()
    for(i in 1:length(cov_data[["obs_t_data"]])){
      times2 = unique(c(times2,cov_data[["obs_t_data"]][[i]][["times"]]))
    }
    times = sort(unique(c(times,times2)))
  }

  covs_movement = unify_movement_covariates(triangulation = domain[["triangulation"]], cov_data = cov_data,
                                            observation_effort = observation_effort, times = times)

  dynamic_covs_movement = (class(covs_movement) == "list")

  if(is.null(cov_data[["obs_t_data"]])){
    dynamic_covs_observation = FALSE
    covs_observation = cov_data[["obs_data"]]
  }else{
    dynamic_covs_observation = TRUE
    covs_observation = cov_data[["obs_t_data"]]
    if(!is.null(cov_data[["obs_data"]])){
      addnames = colnames(cov_data[["obs_data"]])
      for(i in 1:length(covs_observation)){
        x = covs_observation[[i]][["covs"]]
        vals = cov_data[["obs_data"]][i, ]
        for(j in 1:length(addnames)) x[addnames[j]] = vals[j]
        covs_observation[[i]][["covs"]] = x
      }
    }
  }

  site_nodes = matrix(0, ncol = length(observation_effort[["locations"]]), nrow = nrow(domain[["triangulation"]][["node"]]))

  for(i in 1:length(observation_effort[["locations"]])){
    nod_ele = matrix(domain[["triangulation"]][["ele"]][observation_effort[["locations"]][[i]], ], ncol = 3)
    for(j in 1:nrow(nod_ele)){
      site_nodes[nod_ele[j, ], i] = 1
    }
  }

  f_data = affine_mapping(triangulation = domain[["triangulation"]])

  included = rep(NA, 6)
  names(included) = c("diffusion", "advection_1", "advection_2", "mortality", "habitat_preference", "observation")

  for(i in 1:length(included)){
    included[i] = !is.null(model_formula[[names(included)[i]]])
  }

  nc = rep(0, 6)
  names(nc) = names(included)

  X = list()

  for(component in c("diffusion", "advection_1", "advection_2", "mortality", "habitat_preference")){
    if(included[component]){
      if(dynamic_covs_movement){
        tmp = list()
        for(i in 1:length(covs_movement)){
          a = stats::model.matrix(model_formula[[component]], covs_movement[[i]])

          if(component == "habitat_preference"){
            nch = ncol(a)
            b = which(colnames(a) == "(Intercept)")
            if(length(b) == 1){
              name_covs = colnames(a)
              a =  matrix(a[, - b], ncol = (nch - 1))
              colnames(a) = name_covs[-1]
            }
          }
          tmp[[i]] = a
        }
        nc[component] = ncol(tmp[[1]])
      }else{
        a = stats::model.matrix(model_formula[[component]], covs_movement)
        if(component == "habitat_preference"){
          b = which(colnames(a) == "(Intercept)")
          if(length(b) == 1) a = as.matrix(a[, -b], ncol = (ncol(a) - 1) )
        }
        tmp = a
        nc[component] = ncol(tmp)
      }
      X[[component]] = tmp
    }
  }

  observation_effort[["releases"]][["time"]] = match(observation_effort[["releases"]][["time"]], times)
  observation_times = observation_effort[["captures"]][["time"]]

  if(observation_effort[["method"]] == "ICP"){
    observation_times = match(observation_times, times)
  }else{
    for(j in 1:2) observation_times[, j] = match(observation_times[, j], times)
  }
  observation_effort[["captures"]][["time"]] = observation_times

  if(included["observation"]){
    if(dynamic_covs_observation){
      observation_times = list()
      tmp = list()
      for(i in 1:length(covs_observation)){
        tmp[[i]] = stats::model.matrix(model_formula[["observation"]], covs_observation[[i]][["covs"]])
        observation_times[[i]] = match(covs_observation[[i]][["times"]], times)
      }
      nc["observation"] = ncol(tmp[[1]])
    }else{
      if(is.null(covs_observation)){
        tmp = matrix(1, nrow = length(observation_effort[["captures"]][["location"]]))
        colnames(tmp) = "(Intercept)"
      }else{
        tmp = stats::model.matrix(model_formula[["observation"]], covs_observation)
      }
      nc["observation"] = ncol(tmp)
    }
    X[["observation"]] = tmp
    X[["observation_times"]] = observation_times
  }
  if(!(is.null(model_formula[["traits"]]) || is.null(tr_data))){
    X[["trait"]] =  stats::model.matrix(model_formula[["traits"]], tr_data)
  }else{
    tmp = matrix(1, ncol = 1, nrow = ns)
    colnames(tmp) = "(Intercept)"
    X[["trait"]] = tmp
    model_formula[["traits"]] = NULL
  }

  if(is.null(sp_names)){
    sp_names = paste0("species_", 1:ns)
  }

  m = list()
  m[["domain"]] = domain
  m[["observation_effort"]] = observation_effort
  m[["releases"]] = releases
  m[["captures"]] = captures
  m[["cov_data"]] = cov_data
  m[["ns"]] = ns
  m[["sp_names"]] = sp_names
  m[["C"]] = C
  m[["tr_data"]] = tr_data
  m[["model_formula"]] = model_formula
  m[["times"]] = times
  m[["nc"]] = nc
  m[["X"]] = X
  m[["dynamic_covs_observation"]] = dynamic_covs_observation
  m[["dynamic_covs_movement"]] = dynamic_covs_movement
  m[["site_nodes"]] = site_nodes
  m[["covs_movement"]] = covs_movement
  m[["covs_observation"]] = covs_observation
  m[["parNames"]] = labeling_pars(m)
  m[["f_data"]] = f_data
  m[["nt"]] = ncol(m[["X"]][["trait"]])

  de.time = max(m[["times"]]) - min(m[["times"]])
  mu.mortality = log(-log(0.5)/de.time)
  de.space = sqrt(sum(sf::st_area(m[["domain"]][["polygon"]][["geometry"]])))
  e.displacement = de.space / 2
  mu.diffusion = log(exp(mu.mortality) * (2*e.displacement/pi)^2)

  if(m[["observation_effort"]][["method"]] == "CCP"){
    c.times = m[["observation_effort"]][["captures"]][["time"]]
    de.capture.time = mean(m[["times"]][c.times[, 2]] - m[["times"]][c.times[, 1]])
    mu.observation = log(-log(0.5)/(de.capture.time))
  }else{
    mu.observation = 0
  }
  sd.advection = e.displacement*exp(mu.mortality)
  nc = m[["nc"]]
  np = sum(nc)
  nt = m[["nt"]]
  muZ = matrix(0, nrow = nt*np)
  idx = 1
  PSI = diag(np)

  for(j in 1:nt){
    for(i in 1:length(nc)){
      if(names(nc)[i] == "advection_1" && nc[i] > 0){
        for(i_a in 1:nc[i]){
          PSI[(i_a+nc[1]), (i_a + nc[1])] = sd.advection^2
        }
      }
      if(names(nc)[i] == "advection_2" && nc[i] > 0){
        for(i_a in 1:nc[i]){
          PSI[(i_a + (nc[1] + nc[2])), (i_a + (nc[1] + nc[2]))] = sd.advection^2
        }
      }
      #
      for(ii in seq_len(nc[i])){
        muZ[idx] = 0
        if(names(nc)[i] == "diffusion" & ii == 1 & j == 1) muZ[idx] = mu.diffusion
        if(names(nc)[i] == "mortality" & ii == 1 & j == 1) muZ[idx] = mu.mortality
        if(names(nc)[i] == "observation" & ii == 1 & j == 1) muZ[idx] = mu.observation
        idx = idx + 1
      }
    }
  }
  m = set_prior_parameters(m, muZ = muZ, PSI = PSI)
  pre = precomputation_KM(m)

  m[["pre"]] = pre

  class(m) = "class_m"
  return(m)
}
