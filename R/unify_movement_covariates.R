
# R/unify_movement_covariates.R

#' @noRd
#'

unify_movement_covariates = function(triangulation, cov_data, observation_effort, times){
  ne = nrow(triangulation$ele)
  nnod = nrow(triangulation$node)

  is_cov = c(!is.null(cov_data$n_data), !is.null(cov_data$e_data), !is.null(cov_data$t_data),
            !is.null(cov_data$n_t_data), !is.null(cov_data$e_t_data))

  if(any(is_cov[3:5])){
    dt =  times
    ndt = length(dt)
  }

  label  = c()
  label_f = c()

  if(any(is_cov[3:5])){
    D  = list()
    D_f = list()

    for(t in 1:ndt){
      D[[t]] = data.frame(matrix(ncol = 0, nrow = ne*3))
      D_f[[t]] = data.frame(matrix(ncol = 0, nrow = ne*3))
    }
  }else{
    D = data.frame(matrix(ncol = 0, nrow = ne*3))
    D_f = data.frame(matrix(ncol = 0, nrow = ne*3))
  }

  if(is_cov[1]){
    n_id = which(unlist(lapply(cov_data$n_data, is.numeric)) == TRUE)
    f_id = which(unlist(lapply(cov_data$n_data, is.numeric)) == FALSE)

    if(length(n_id) > 0){
      D_aux = spatial_ext(triangulation = triangulation, s_data = data.frame(matrix(cov_data$n_data[, n_id])) )
    }
    if(length(f_id) > 0){
      D_aux_f = spatial_ext_f(triangulation = triangulation, s_data = data.frame(matrix(cov_data$n_data[, f_id])) )
    }

    if(any(is_cov[3:5])){
      for(t in 1:ndt){
        if(length(n_id) > 0){
          D[[t]] = D_aux
        }
        if(length(f_id) > 0){
          D_f[[t]] = D_aux_f
        }
      }
    }else{
      if(length(n_id) > 0){
        D = D_aux
      }
      if(length(f_id) > 0){
        D_f = D_aux_f
      }
    }
    label  = names(n_id)
    label_f = names(f_id)
  }

  if(is_cov[2]){
    n_id = which(unlist(lapply(cov_data$e_data, is.numeric)) == TRUE)
    f_id = which(unlist(lapply(cov_data$e_data, is.numeric)) == FALSE)

    if(length(n_id) > 0){
      D_aux = data.frame(matrix(cov_data$e_data[rep(1:ne, 3), n_id]))
    }

    if(length(f_id) > 0){
      D_aux_f = data.frame(matrix(cov_data$e_data[rep(1:ne, 3), f_id]))
    }

    if(any(is_cov[3:5])){
      for(t in 1:ndt){
        if(length(n_id) > 0){
          D[[t]] = cbind(D[[t]], D_aux)
        }
        if(length(f_id) > 0){
          D_f[[t]] = cbind(D_f[[t]], D_aux_f)
        }
      }
    }else{
      if(length(n_id) > 0){
        D = cbind(D, D_aux)
      }
      if(length(f_id) > 0){
        D_f = cbind(D_f, D_aux_f)
      }
    }

    label  = c(label, names(n_id))
    label_f = c(label_f, names(f_id))
  }

  if(is_cov[3]){
    n_id = which(unlist(lapply(cov_data$t_data, is.numeric)) == TRUE)
    f_id = which(unlist(lapply(cov_data$t_data, is.numeric)) == FALSE)

    times_l  = cov_data$t_data_times
    s_t_data  = list()
    s_t_data_f = list()

    for(t in 1:length(times_l)){
      if(length(n_id) > 0){
        if(length(n_id)==1){
          s_t_data[[t]] = data.frame(matrix(cov_data$t_data[rep(t, ne*3), n_id], ncol = length(n_id) ))
        }else{
          s_t_data[[t]] = cov_data$t_data[rep(t, ne*3), n_id]
        }
      }
      if(length(f_id) > 0){
        s_t_data_f[[t]] = data.frame(cov_data$t_data[rep(t, ne*3), f_id])
      }
    }

    D_aux  = interpola_t(s_t_data = s_t_data, s_t_data_times = times_l, dt)
    D_aux_f = interpola_t_f(s_t_data = s_t_data_f, s_t_data_times = times_l, dt)

    for(t in 1:ndt){
      if(length(n_id) > 0){
        D[[t]] = cbind(D[[t]], D_aux[[t]])
      }
      if(length(f_id) > 0){
        D_f[[t]] = cbind(D_f[[t]], D_aux_f[[t]])
      }
    }

    label  = c(label, names(n_id))
    label_f = c(label_f, names(f_id))
  }

  if(is_cov[4]){
    n_id = which(unlist(lapply(cov_data$n_t_data[[1]], is.numeric)) == TRUE)
    f_id = which(unlist(lapply(cov_data$n_t_data[[1]], is.numeric)) == FALSE)

    times_l = cov_data$n_t_data_times
    s_t_data  = list()
    s_t_data_f = list()

    for(t in 1:length(times_l)){
      if(length(n_id) > 0){
        s_t_data[[t]] = spatial_ext(triangulation = triangulation,
                                     s_data = matrix(cov_data$n_t_data[[t]][, n_id], ncol = length(n_id)))
      }
      if(length(f_id) > 0){
        s_t_data_f[[t]] = spatial_ext_f(triangulation = triangulation,
                                        s_data = matrix(cov_data$n_t_data[[t]][, f_id], ncol = length(f_id)))
      }
    }

    D_aux  = interpola_t(s_t_data = s_t_data, s_t_data_times = times_l, dt)
    D_aux_f = interpola_t_f(s_t_data = s_t_data_f, s_t_data_times = times_l, dt)

    for(t in 1:ndt){
      if(length(n_id) > 0){
        D[[t]] = cbind(D[[t]], D_aux[[t]])
      }
      if(length(f_id) > 0){
        D_f[[t]] = cbind(D_f[[t]], D_aux_f[[t]])
      }
    }

    label  = c(label, names(n_id))
    label_f = c(label_f, names(f_id))
  }

  if(is_cov[5]){
    n_id = which(unlist(lapply(cov_data$e_t_data[[1]], is.numeric)) == TRUE)
    f_id = which(unlist(lapply(cov_data$e_t_data[[1]], is.numeric)) == FALSE)

    times_l = cov_data$e_t_data_times
    s_t_data = list()
    s_t_data_f = list()

    for(t in 1:length(times_l)){
      if(length(n_id) > 0){
        if(length(n_id)==1){
          s_t_data[[t]] = data.frame(matrix(cov_data$e_t_data[[t]][rep(1:ne, 3), n_id], ncol = length(n_id)))
        }else{
          s_t_data[[t]] = cov_data$e_t_data[[t]][rep(1:ne, 3), n_id]
        }
      }
      if(length(f_id) > 0){
        if(length(f_id)==1){
          s_t_data_f[[t]] = data.frame(matrix(cov_data$e_t_data[[t]][rep(1:ne, 3), f_id], ncol = length(f_id)))
        }else{
          s_t_data_f[[t]] = cov_data$e_t_data[[t]][rep(1:ne, 3), f_id]
        }
      }
    }

    D_aux  = interpola_t(s_t_data = s_t_data, s_t_data_times = times_l, dt)
    D_aux_f = interpola_t_f(s_t_data = s_t_data_f, s_t_data_times = times_l, dt)

    for(t in 1:ndt){
      if(length(n_id) > 0){
        D[[t]] = cbind(D[[t]], D_aux[[t]])
      }
      if(length(f_id) > 0){
        D_f[[t]] = cbind(D_f[[t]], D_aux_f[[t]])
      }
    }

    label  = c(label, names(n_id))
    label_f = c(label_f, names(f_id))
  }

  if(length(label) > 0 && length(label_f) > 0){
    if(any(is_cov[3:5])){
      for(t in 0:(length(dt) - 1)){
        D[[t + 1]] = cbind(D[[t + 1]], D_f[[t + 1]])

        colnames(D[[t + 1]]) = c(label, label_f)
      }
    }else{
      D = cbind(D,D_f)

      colnames(D) = c(label, label_f)
    }
  }

  if(length(label) > 0 && length(label_f) == 0 ){
    if(any(is_cov[3:5])){
      for(t in 0:(length(dt) - 1)){
        colnames(D[[t + 1]]) = label
      }
    }else{
      colnames(D) = label
    }
  }

  if(length(label) == 0 && length(label_f) > 0 ){
    D = D_f
    if(any(is_cov[3:5])){
      for(t in 0:(length(dt) - 1)){
        colnames(D[[t + 1]]) = label_f
      }
    }else{
      colnames(D) = label_f
    }
  }

  return(D)
}
