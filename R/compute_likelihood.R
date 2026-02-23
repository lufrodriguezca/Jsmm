# R/compute_likelihood.R

#' Compute the likelihood
#'
#' Compute the likelihood in log form for the Joint species movement model
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (jsmm) analyses.
#' @param Beta A named list with the movement and observation parameters included in the Jsmm model.
#' @param max_dt Optional positive integer value representing the number of time partitions in which the PDE should be solved.
#'
#' @return A named list with the loglikelihood computations for each species.
#' @export
#'
#'
compute_likelihood = function(m, Beta, max_dt = NULL){
  try.error = rep(FALSE, m[["ns"]])
  pre = m[["pre"]]
  vectorized = TRUE
  eps = 1e-8
  site_nodes = list()

  for(i in 1:ncol(m[["site_nodes"]])){
    site_nodes[[i]] = which(m[["site_nodes"]][,i] == 1)
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
  numerical_problems = 0
  numerical_attempts = 0
  releases = m[["merged.releases"]]
  captures = m[["merged.captures"]]
  ni = nrow(releases)
  nnod = nrow(m[["domain"]][["triangulation"]][["node"]])
  v = matrix(0, nrow = nnod, ncol = ni)
  release_times = m[["observation_effort"]][["releases"]][["time"]][releases[, 2]]
  release_sites = m[["observation_effort"]][["releases"]][["location"]][releases[, 2]]
  release_species = releases[, 1]
  released = rep(FALSE,ni)
  log.likelihood = rep(0, ni)

  re.dt.mode = -1
  if(!time_dependent_parameters){
    XBeta = multiply_XBeta(m = m, Beta = Beta, 1)
    K = compute_Stiffness(m = m, XBeta = XBeta, pre)
    M = compute_Mass(m = m, XBeta = XBeta, pre)
    re.dts = c()
    for(i1 in 1:(length(m[["times"]]) - 1)){
      time1 = m[["times"]][i1]
      time2 = m[["times"]][i1 + 1]
      ltimes = seq(from = time1,to = time2,length.out = 1 + ceiling((time2 - time1)/max_dt))

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

  if(m[["observation_effort"]][["method"]] == "CCP"){
    p = matrix(nrow = length(m[["observation_effort"]][["captures"]][["location"]]), ncol = ni)
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
        v[, ind] = 0
        nodes = site_nodes[[release_sites[ind]]]
        v[nodes, ind] = 1/sum(integral[nodes, release_species[ind]])
        released[ind] = TRUE
      }
    }

    time1 = m[["times"]][i1]
    time2 = m[["times"]][i1 + 1]
    ltimes = seq(from = time1,to = time2,length.out = 1 + ceiling((time2 - time1)/max_dt))

    if(m[["observation_effort"]][["method"]] == "CCP"){
      p[which(m[["observation_effort"]][["captures"]][["time"]][, 1] == i1), ] = 0
      capture_events = which(m[["observation_effort"]][["captures"]][["time"]][, 1]<= (i1) & m[["observation_effort"]][["captures"]][["time"]][, 2] >= (i1 + 1))
    }

    for(i2 in 1:(length(ltimes)-1)){
      current_time = (ltimes[i2] + ltimes[i2+1])/2
      re.dt = ltimes[i2 + 1] - ltimes[i2]
      if(time_dependent_parameters){
        w2 = (current_time-time1)/(time2 - time1)
        w1 = 1 - w2
        XBeta = lapply(XBeta_start, "*", w1)
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
            D2.x = M.x  + (1-(2/3))*re.dt*K.x
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
      inds =  which(released)
      spp = release_species[inds]
      for(sp in unique(spp)){
        if(!try.error[sp]){
          ind = inds[release_species[inds] == sp]
          v_old = matrix(v[,ind], ncol = length(ind))
          u = D2[[sp]]%*%v[, ind]
          v_new = try(matrix(Matrix::solve(D1[[sp]], u), ncol = length(ind)), silent = TRUE)
          if(any(class(v_new) == "try-error")){
            try.error[sp] = TRUE
            v_new = v_old
          }
          if(any(is.na(v_new))){
            try.error[sp] = TRUE
            v_new = v_old
          }
          if(any(is.infinite(v_new))){
            try.error[sp] = TRUE
            v_new = v_old
          }
          v[,ind] = v_new
          if(m[["observation_effort"]][["method"]] == "CCP" && length(capture_events) > 0){
            if(m[["dynamic_covs_observation"]]){
              r = matrix(0,nrow = length(capture_events), ncol = m[["ns"]])
              for(cc1 in 1:length(capture_events)){
                cc = capture_events[cc1]
                nn = which(current_time < m[["times"]][m[["X"]][["observation_times"]][[cc]]])[1] - 1
                r[cc1, ] = m[["X"]][["observation"]][[cc]][nn, ]%*% t(Beta[["observation"]])
              }
            } else {
              r = matrix(m[["X"]]$observation[capture_events, ], nrow = length(capture_events)) %*% t(Beta[["observation"]])
            }
            r = exp(r)
            v_mean = (v_old+v[, ind])/2
            for(cc1 in 1:length(capture_events)){
              cc = capture_events[cc1]
              capture_site = m[["observation_effort"]][["captures"]][["location"]][cc]
              a1 = integral[site_nodes[[capture_site]], sp]
              a2 = matrix(v_mean[site_nodes[[capture_site]], ],ncol = length(ind))
              p[cc, ind] = p[cc, ind] + re.dt*colSums(a1*a2)*r[cc1, sp]
            }
          }
        }
      }
    }
    if(time_dependent_parameters){
      XBeta = XBeta_end
      if(time_dependent_k){
        integral = compute_integral(m = m, XBeta = XBeta)
      }
    }

    if(m[["observation_effort"]][["method"]] == "CCP"){
      capture_closes = which(m[["observation_effort"]][["captures"]][["time"]][, 2] == i1 + 1)
      all.captures.cc = list()
      for(cc1 in seq_len(length(capture_closes))){
        all.captures.cc[[cc1]] = matrix(captures[(captures[, 2] == capture_closes[cc1]), ], ncol = 3)
      }
      if(length(capture_closes) > 0){
        for(ind in which(released)){
          captured_this_time = FALSE
          for(cc1 in seq_len(length(capture_closes))){
            cc = capture_closes[cc1]
            captures.cc = all.captures.cc[[cc1]]
            capture_site = m[["observation_effort"]][["captures"]][["location"]][cc]
            if(!captured_this_time){
              qq = min(1 - eps, max(eps, p[cc, ind]))
              fc = which((captures.cc[, 1] == ind))
              numerical_attempts = numerical_attempts + 1
              if(length(fc) == 0){
                log.likelihood[ind] = log.likelihood[ind] + log(1 - qq)
                v[, ind] = v[, ind]/(1 - qq)
              } else {
                log.likelihood[ind] = log.likelihood[ind] + log(qq)
                captured_this_time = TRUE
                secondary_release = captures.cc[fc, 3]
                if(secondary_release){
                  v[, ind] = 0
                  nodes = site_nodes[[capture_site]]
                  v[nodes, ind] = 1/sum(integral[nodes, release_species[ind]])
                  released[ind] = TRUE
                } else{
                  released[ind] = FALSE
                }
              }
            }
          }
        }
      }
    }


    if(m[["observation_effort"]][["method"]] == "ICP"){
      capture_events = which(m[["observation_effort"]][["captures"]][["time"]] == (i1 + 1))
      if(length(capture_events) > 0){
        p = matrix(m[["X"]]$observation[capture_events, ], nrow = length(capture_events)) %*% t(Beta[["observation"]])
        p = exp(p)/(1 + exp(p))
        for(cc1 in 1:length(capture_events)){
          cc = capture_events[cc1]
          capture_site = m[["observation_effort"]][["captures"]][["location"]][cc]
          captures.cc = matrix(captures[(captures[, 2] == cc), ], ncol = 3)
          for(ind in which(released)){
            sp = release_species[ind]
            U = sum(integral[site_nodes[[capture_site]], sp]*v[site_nodes[[capture_site]], ind])
            qq = min(1 - eps,max(eps,p[cc1,sp]*U))
            numerical_attempts = numerical_attempts + 1
            fc = which((captures.cc[, 1] == ind))
            if(length(fc) == 0){
              if(qq >= 1){
                log.likelihood[ind] = log.likelihood[ind] + log(eps)
                numerical_problems = numerical_problems + 1
              } else {
                log.likelihood[ind] = log.likelihood[ind] + log(1 - qq)
                v[, ind] = v[, ind]/(1 - qq)
              }
              nodes = site_nodes[[capture_site]]
              v[nodes, ind] = v[nodes, ind]*(1 - p[cc1, sp])
            } else {
              if(qq <= 0){
                log.likelihood[ind] = log.likelihood[ind] + log(eps)
                numerical_problems = numerical_problems + 1
              } else{
                log.likelihood[ind] = log.likelihood[ind] + log(qq)
              }
              secondary_release = captures.cc[fc, 3]
              if(secondary_release){
                v[, ind] = 0
                nodes = site_nodes[[capture_site]]
                v[nodes, ind] = 1/sum(integral[nodes, release_species[ind]])
              } else{
                released[ind] = FALSE
              }
            }
          }
        }
      }
    }
  }

  log.likelihood.sp = rep(NA,m[["ns"]])
  for(i in 1:m[["ns"]]){
    sel = which(release_species == i)
    log.likelihood.sp[i] = sum(log.likelihood[sel]*releases[sel, 3])
    if(try.error[i])   log.likelihood.sp[i] = - 10^6
  }
  res = list()
  res[["log.likelihood.sp"]] = log.likelihood.sp
  res[["numerical_problems"]] = numerical_problems
  res[["try.error"]] = try.error
  return(res)
}

