# R/spatial_ext.R

#' @noRd
#'
spatial_ext = function(triangulation, s_data){
  nele = nrow(triangulation$ele)
  s_data_ext = data.frame(matrix(NA, ncol = ncol(s_data), nrow = nele*3))

  for(i in 1:nele){
    s_data_ext[i, ] = s_data[triangulation$ele[i, 1], ]
    s_data_ext[i + nele, ] = s_data[triangulation$ele[i, 2], ]
    s_data_ext[i + (nele*2), ] = s_data[triangulation$ele[i, 3], ]
  }

  s_data_ext = data.frame(s_data_ext)

  return(s_data_ext)
}

