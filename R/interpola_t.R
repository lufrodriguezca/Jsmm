# R/interpola_t.R

#' @noRd
#'
interpola_t = function(s_t_data, s_t_data_times, dt){
  s_t_data_ext = list()

  if(length(s_t_data) > 0){
    nt  = length(s_t_data_times)
    ncov = ncol(s_t_data)

    for(t in 1:(nt - 1)){
      y1 = s_t_data[[t]]
      y2 = s_t_data[[t + 1]]
      x1 = s_t_data_times[[t]]
      x2 = s_t_data_times[[t + 1]]
      m = (y2 - y1)/(x2 - x1)
      b = y1 - m*x1

      subinterval = which((dt > s_t_data_times[t] | dt == s_t_data_times[t]) &
                          (dt < s_t_data_times[t + 1] | dt == s_t_data_times[t + 1]))

      for(k in subinterval){
        s_t_data_ext[[k]] = m*dt[k] + b
      }
    }
  }

  return(s_t_data_ext)
}
