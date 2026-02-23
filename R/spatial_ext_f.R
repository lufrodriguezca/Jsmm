# R/spatial_ext_f.R

#' @noRd
#'
spatial_ext_f = function(triangulation, s_data){
  nele = nrow(triangulation$ele)
  s_data_ext = data.frame(matrix(NA, ncol = ncol(s_data), nrow = nele*3))

  for(i in 1:nele){
    s_data_ext[i, ] = as.character(s_data[triangulation$ele[i, 1], ])
    s_data_ext[i + (nele), ] = as.character(s_data[triangulation$ele[i, 2], ])
    s_data_ext[i + (nele*2), ] = as.character(s_data[triangulation$ele[i, 3], ])
  }

  s_data_ext = data.frame(lapply(s_data_ext, factor))

  return(s_data_ext)
}

