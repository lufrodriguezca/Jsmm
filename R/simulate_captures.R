# R/simulate_captures.R

#' Simulate movement data.
#'
#' The function returns the model object m, to which it has added the objects captures and releases.
#'
#' @param m named list created by the function jsmm().
#' @param pars named list with the objects Beta, rho, SI, and Gamma defining the model parameters.
#' @param releases matrix whose first and second columns provide information about the species and release event, respectively.
#' @param max_dt real number bigger than zero that describes the maximal time resolution used in numerical computations.
#' @param secondary_release boolean parameter that describes whether the captured individuals are released after capture.
#'
#' @return The function returns the model object m, to which it has added the objects captures and releases.
#' @export
#'

simulate_captures = function(m, pars, releases,max_dt = NULL, secondary_release){
  vectorized = TRUE
  Beta = pars[["Beta"]]
  site_nodes = list()
  for(i in 1:ncol(m[["site_nodes"]])){
    site_nodes[[i]] =   which(m[["site_nodes"]][, i] == 1)
  }
  time_dependent_parameters = (m[["dynamic_covs_observation"]] ||  m[["observation_effort"]][["method"]]=="CCP" ||  m[["dynamic_covs_movement"]])
  time_dependent_k = FALSE
  if(time_dependent_parameters){
    kvars = all.vars(m[["model_formula"]][["habitat_preference"]])
    dynvars = colnames(m[["cov_data"]][["t_data"]])
    if(!is.null(m[["cov_data"]][["e_t_data"]])){
      dynvars = c(dynvars,colnames(m[["cov_data"]][["e_t_data"]][[1]]))
    }
    time_dependent_k = any(kvars %in% dynvars)
  }
  if(is.null(max_dt)) max_dt = (max(m[["times"]]) - min(m[["times"]]))/50

  pre = m[["pre"]]
  ni = nrow(releases)

  captures = c()
  nnod = nrow(m[["domain"]][["triangulation"]][["node"]])
  v = matrix(0, nrow = nnod, ncol = ni)
  release_times = m[["observation_effort"]][["releases"]][["time"]][releases[, 2]]
  release_sites = m[["observation_effort"]][["releases"]][["location"]][releases[, 2]]
  release_species = releases[, 1]
  released_and_moving = rep(FALSE, ni)

  re.dt.mode = -1
  if(!time_dependent_parameters){
    XBeta = multiply_XBeta(m = m, Beta = Beta, 1)
    K = compute_Stiffness(m = m, XBeta = XBeta, pre)
    M = compute_Mass(m = m, XBeta = XBeta, pre)
    re.dts = c()
    for(i1 in 1:(length(m[["times"]]) - 1)){
      time1 = m[["times"]][i1]
      time2 = m[["times"]][i1 + 1]
      ltimes = seq(from=time1,to = time2,length.out = 1 + ceiling((time2 - time1)/max_dt))
      for(i2 in 1:(length(ltimes) - 1)){
        re.dt = ltimes[i2 + 1] - ltimes[i2]
        re.dts = c(re.dts, re.dt)
      }
    }
    re.dt.mode = as.numeric(names(sort(table(re.dts), decreasing = T))[1])
    D1.mode = list()
    D2.mode = list()
    for(sp in 1:m[["ns"]]){
      if(vectorized){
        M.x = t(M[[sp]])@x
        K.x = t(K[[sp]])@x
        D1.x = M.x - (2/3)*re.dt.mode*K.x
        D2.x = M.x  + (1 - (2/3))*re.dt.mode*K.x
        D1.mode[[sp]] = t(M[[sp]])
        D2.mode[[sp]] = t(M[[sp]])
        D1.mode[[sp]]@x = D1.x
        D2.mode[[sp]]@x = D2.x}
      else {
        D1.mode[[sp]] = t(M[[sp]]) - (2/3)*re.dt.mode*t(K[[sp]])
        D2.mode[[sp]] = t(M[[sp]]) + (1-(2/3))*re.dt.mode*t(K[[sp]])
      }
    }
  }
  if(!time_dependent_k){
    XBeta = multiply_XBeta(m = m, Beta = Beta, 1)
    integral = compute_integral(m = m, XBeta = XBeta)
  }

  for(i1 in 1:(length(m[["times"]]) - 1)){
    if(time_dependent_parameters){
      if(i1 == 1){
        XBeta_start = multiply_XBeta(m = m, Beta = Beta, i = 1)
        XBeta_end = multiply_XBeta(m = m, Beta = Beta, i = 2)
      } else {
        XBeta_start = XBeta_end
        XBeta_end = multiply_XBeta(m = m, Beta = Beta, i = (i1 + 1))
      }
    }

    released_inds = which(release_times == i1)
    if(length(released_inds) > 0){
      if(time_dependent_parameters){
        XBeta = XBeta_start
      }
      if(time_dependent_k){
        integral = compute_integral(m = m, XBeta = XBeta)
      }
      for(ind in released_inds){
        nodes = site_nodes[[release_sites[ind]]]
        current_node = sample(nodes, prob = integral[nodes, release_species[ind]], size = 1)
        v[, ind] = 0
        v[current_node, ind] = 1/integral[current_node, release_species[ind]]
        released_and_moving[ind] = TRUE
      }
    }

    time1 = m[["times"]][i1]
    time2 = m[["times"]][i1 + 1]
    delta = time2-time1
    ltimes = seq(from = time1, to = time2, length.out = 1+ceiling(delta/max_dt))

    if(m[["observation_effort"]][["method"]] == "CCP"){
      active_capture_events = which(i1>=m[["observation_effort"]][["captures"]][["time"]][,1] & (i1+1)<=m[["observation_effort"]][["captures"]][["time"]][,2]) #see which capture events are active
      p = matrix(0,nrow=length(active_capture_events),ncol = ni)
      if(length(active_capture_events) > 0 & !m[["dynamic_covs_observation"]]){
        r = matrix(m[["X"]][["observation"]][active_capture_events, ], nrow = length(active_capture_events)) %*% t(Beta[["observation"]])
        r = exp(r)
      }
    }

    for(i2 in 1:(length(ltimes) - 1)){
      current_time = (ltimes[i2] + ltimes[i2 + 1])/2
      re.dt = ltimes[i2 + 1] - ltimes[i2]
      if(time_dependent_parameters){
        w2 = (current_time-time1)/(time2 - time1)
        w1 = 1 - w2
        XBeta = lapply(XBeta_start,"*", w1)
        tmp = lapply(XBeta_end, "*", w2)
        for(i in 1:length(XBeta)) XBeta[[i]] = XBeta[[i]] + tmp[[i]]
        K = compute_Stiffness(m = m, XBeta = XBeta, pre)
        M = compute_Mass(m = m, XBeta = XBeta, pre)
      }
      if(time_dependent_k){
        integral = compute_integral(m = m, XBeta = XBeta)
      }
      if(!time_dependent_parameters & re.dt == re.dt.mode){
        D1 = D1.mode
        D2 = D2.mode
      } else {
        D1 = list()
        D2 = list()
        for(sp in 1:m[["ns"]]){
          if(vectorized){
            M.x = t(M[[sp]])@x
            K.x = t(K[[sp]])@x
            D1.x = M.x - (2/3)*re.dt*K.x
            D2.x = M.x  + (1 - (2/3))*re.dt*K.x
            D1[[sp]] = t(M[[sp]])
            D2[[sp]] = t(M[[sp]])
            D1[[sp]]@x = D1.x
            D2[[sp]]@x = D2.x}
          else {
            D1[[sp]] = t(M[[sp]]) - (2/3)*re.dt*t(K[[sp]])
            D2[[sp]] = t(M[[sp]]) + (1 - (2/3))*re.dt*t(K[[sp]])
          }
        }
      }
      if(m[["observation_effort"]][["method"]] == "CCP"){
        if(length(active_capture_events) > 0 & m[["dynamic_covs_observation"]]){
          r = matrix(0, nrow = length(active_capture_events), ncol = m[["ns"]])
          for(cc in 1:length(active_capture_events)){
            nn = which(current_time < m[["times"]][m[["X"]][["observation_times"]][[active_capture_events[cc]]]])[1] - 1
            r[cc, ] = m[["X"]][["observation"]][[cc]][nn, ] %*% t(Beta[["observation"]])
          }
          r = exp(r)
        }
      }

      all.inds =  which(released_and_moving)
      spp = release_species[all.inds]
      v_old = v
      for(sp in unique(spp)){
        inds = all.inds[release_species[all.inds] == sp]
        u = D2[[sp]] %*% v[, inds]
        v[, inds] = matrix(Matrix::solve(D1[[sp]], u), ncol = length(inds))
        if(m[["observation_effort"]][["method"]] == "CCP" && length(active_capture_events) > 0){
          for(ind in inds){
            v_mean = (v_old[, ind] + v[, ind])/2
            for(cc1 in 1:length(active_capture_events)){
              cc = active_capture_events[cc1]
              capture_site = m[["observation_effort"]][["captures"]][["location"]][cc]
              p[cc1, ind] =  p[cc1, ind] + r[cc1, sp]*re.dt*sum(integral[site_nodes[[capture_site]], sp]*v_mean[site_nodes[[capture_site]]])
            }
          }
        }
      }
    }

    if(m[["observation_effort"]][["method"]] == "CCP"){
      for(ind in which(released_and_moving)){
        p_alive_and_not_captured = sum(v[,ind]*integral[, release_species[ind]])
        if(length(active_capture_events) > 0){
          p_captured = sum(p[, ind])
          p_dead_or_captured = 1 - p_alive_and_not_captured
          p_dead = p_dead_or_captured - p_captured
        } else {
          p_dead = 1 - p_alive_and_not_captured
        }
        if(runif(1) < p_dead){
          released_and_moving[ind] = FALSE
        } else {
          v[, ind] = v[, ind]/(1 - p_dead)
          if(length(active_capture_events) > 0){
            p[, ind] = p[, ind]/(1 - p_dead)
            if(runif(1) < sum(p[, ind])){
              ce = active_capture_events[sample(1:length(p[, ind]), prob = pmax(0, p[, ind]), size = 1)]
              captures = rbind(captures,c(ind, ce, secondary_release))
              released_and_moving[ind] = secondary_release
              if(secondary_release){
                nodes = site_nodes[[m[["observation_effort"]][["captures"]][["location"]][ce]]]
                current_node = sample(nodes, prob = integral[nodes, release_species[ind]], size = 1)
                v[, ind] = 0
                v[current_node, ind] = 1/integral[current_node, release_species[ind]]
              }
            } else {
              v[, ind] = v[, ind]/sum(v[, ind]*integral[, release_species[ind]])
            }
          }
        }
      }
    }

    if(m[["observation_effort"]][["method"]] == "ICP"){
      active_capture_events = which(m[["observation_effort"]][["captures"]][["time"]] == (i1 + 1))
      if(length(active_capture_events) > 0){
        p = matrix(m[["X"]]$observation[active_capture_events, ], nrow = length(active_capture_events)) %*% t(Beta[["observation"]])
        p = exp(p)/(1 + exp(p))
        if(time_dependent_k){
          if(time_dependent_parameters){
            integral = compute_integral(m = m, XBeta = XBeta_end)
          } else {
            integral = compute_integral(m = m, XBeta = XBeta)
          }
        }
      }

      for(ind in which(released_and_moving)){
        p_alive =  sum(v[, ind]*integral[, release_species[ind]])
        if(runif(1) > p_alive){
          released_and_moving[ind] = FALSE
        } else {
          v[, ind] = v[, ind]/p_alive
          pp = rep(length(active_capture_events))
          for(cc1 in 1:length(active_capture_events)){
            cc = active_capture_events[cc1]
            capture_site = m[["observation_effort"]][["captures"]][["location"]][cc]
            nodes = site_nodes[[capture_site]]
            pp[cc1] = p[cc1, release_species[ind]]*sum(v[nodes, ind]*integral[nodes, release_species[ind]])
          }
          if(runif(1) < sum(pp)){
            ce = active_capture_events[sample(1:length(pp), prob = pmax(0, pp), size = 1)]
            captures = rbind(captures, c(ind, ce, secondary_release))
            released_and_moving[ind] = secondary_release
            if(secondary_release){
              nodes = site_nodes[[m[["observation_effort"]][["captures"]][["location"]][ce]]]
              current_node = sample(nodes, prob = integral[nodes, release_species[ind]], size = 1)
              v[, ind] = 0
              v[current_node, ind] = 1/integral[current_node, release_species[ind]]
            }
          } else {
            v[, ind] = v[, ind]/sum(v[, ind]*integral[, release_species[ind]])
          }
        }
      }
    }
  }

  if(!is.null(captures)){
    colnames(captures) = c("individual", "capture_event", "secondary_release")
  }

  m1 = m
  m1$releases = releases
  m1$captures = captures

  class(m1) = "class_m"
  return(m1)
}


