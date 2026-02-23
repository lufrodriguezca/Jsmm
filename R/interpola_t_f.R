# R/interpola_t_f.R

#' @noRd
#'
interpola_t_f = function(s_t_data, s_t_data_times, dt){
  s_t_data_ext_f    = list()
  if(length(s_t_data) > 0){
    nt  = length(s_t_data_times)
    s_t_data_ext_f = list()

    for(t in 1:(nt - 1)){
      subinterval_f = which((dt > s_t_data_times[t] | dt == s_t_data_times[t]) &
                            (dt < s_t_data_times[t + 1] | dt == s_t_data_times[t + 1]))

      for(k in subinterval_f){
        e1 = abs(s_t_data_times[t] - dt[k] )
        e2 = abs(s_t_data_times[t + 1] - dt[k])

        if(e1 < e2){
          s_t_data_ext_f[[k]] = s_t_data[[t]]
        }
        if(e1 > e2){
          s_t_data_ext_f[[k]] = s_t_data[[t + 1]]
        }
        if(e1 == e2){
          s_t_data_ext_f[[k]] = s_t_data[[sample(c(t, t + 1), 1)]]
        }
      }
    }
  }

  return(s_t_data_ext_f)
}

