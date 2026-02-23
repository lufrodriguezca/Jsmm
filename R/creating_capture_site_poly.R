# R/creating_capture_site_poly.R

#' @noRd
#'

creating_capture_site_poly = function(n, c, rad, theta){
  coor_t = matrix(NA, ncol = 2, nrow = (n + 1))

  for(i in 1:n){
	  coor_t[i, 1] = rad*cos(2*pi*i/n + theta) + c[1]
		coor_t[i, 2] = rad*sin(2*pi*i/n + theta) + c[2]
  }

  coor_t[(n + 1), ] = coor_t[1, ]

  return(coor_t)
}
