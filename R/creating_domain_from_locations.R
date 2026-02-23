#source("functions/creating_capture_site_poly.R")

# R/creating_domain_from_locations.R

#' Create a simulated domain
#'
#' This function allows the creation of a customized rectangular domain, capture, and release sites.
#'
#' @param lc matrix that provides the positions of the capture and release sites. Columns correspond to the x and y coordinates, respectively.
#' @param lc_label character vector with capture/release sites labels (optional).
#' @param buffer_x real positive number representing the distance from the center of the outer capture/release location to the domain boundary in the x-axis.
#' @param buffer_y real positive number representing the distance from the center of the outer capture/release location to the domain boundary in the y-axis.
#' @param  nnod_side_x integer positive representing the number of nodes in x axes.
#' @param  nnod_side_y integer positive representing the number of nodes in y axes.
#' @param  rad_capture_site real positive value representing the capture/release site radius.
#' @param  n_sides_capture_site Number of sides for the capture/release polygon.
#' @param  angle_capture_site input parameter in radians measuring the counterclockwise rotation angle for the capture/release polygon.
#' @param  int_buffer boolean parameter for inserting an intermediate buffer between the border of the captures/releases array and the landscape boundary.
#'
#' @return Polygonal object from the class sf.
#' @export
#'
creating_domain_from_locations = function(lc, lc_label = NA, buffer_x = c(0,0) , buffer_y = c(0, 0), nnod_side_x,nnod_side_y,
                                   rad_capture_site, n_sides_capture_site, angle_capture_site, int_buffer = FALSE){

  if(!exists(deparse(substitute(lc)))){
    message("STOP: Please insert coordinates of locations.")
  }else{
    if(sum(is.na(lc_label))>0){
      lc_label = 1:nrow(lc)
    }
  }

  minx = min(lc[, 1])
  miny = min(lc[, 2])
  maxx = max(lc[, 1])
  maxy = max(lc[, 2])

  xnod = seq((minx - buffer_x[1]), (maxx + buffer_x[2]), by = (((maxx + buffer_x[2]) - (minx - buffer_x[1]))/(nnod_side_x - 1)))
  ynod = seq((miny - buffer_y[1]), (maxy + buffer_y[2]), by = (((maxy + buffer_y[2]) - (miny - buffer_y[1]))/(nnod_side_y - 1)))

  edge1 = matrix(NA, ncol = 2, nrow = (nnod_side_y - 1) )
  edge1[,1] = (minx - buffer_x[1])
  edge1[,2] = sort(ynod,decreasing=TRUE)[-nnod_side_y]

  edge2 = matrix(NA, ncol = 2, nrow = (nnod_side_x - 1)  )
  edge2[,1] = xnod[-nnod_side_x]
  edge2[,2] = (miny-buffer_y[1])

  edge3 = matrix(NA, ncol = 2, nrow = (nnod_side_y-1)  )
  edge3[,1] = (maxx + buffer_x[2])
  edge3[,2] = ynod[-nnod_side_y]

  edge4 = matrix(NA, ncol = 2, nrow = nnod_side_x )
  edge4[,1] = sort(xnod,decreasing=TRUE)
  edge4[,2] = (maxy + buffer_y[2])
  polycoor = rbind(edge1, edge2, edge3, edge4)

  border_dom = list()
  border_dom[[1]] = st_polygon(list(polycoor))
  border_dom = st_sf(id = 1, id_location = 0, geometry = st_sfc(border_dom))

  border_dom_repeat = list()
  border_dom_repeat[[1]] = st_polygon(list(polycoor))
  border_dom_repeat = st_sf(id = 2, id_location = 0, geometry = st_sfc(border_dom_repeat))

  sites_l = vector("list", (nrow(lc)))

  for(i in 1:nrow(lc)){
    poly_trap = creating_capture_site_poly(n = n_sides_capture_site, c = lc[i,],
                                           rad = rad_capture_site, theta = angle_capture_site)

    sites_l[[i]] = st_polygon(list(poly_trap))
  }

  s = st_sf(id = rep(2, nrow(lc)), id_location = 1:nrow(lc),  geometry = st_sfc(sites_l))

  if( int_buffer && ((minx - buffer_x[1]/2) > min(st_coordinates( s)[,1]) || (minx - buffer_x[1]/2) == min(st_coordinates( s)[, 1]) ||
                     (maxx + buffer_x[2]/2) < max(st_coordinates( s)[,1]) || (maxx + buffer_x[2]/2) == max(st_coordinates( s)[, 1]) ||
                     (miny - buffer_y[1]/2) > min(st_coordinates( s)[,2]) || (miny - buffer_y[1]/2) == min(st_coordinates( s)[, 2]) ||
                     (maxy + buffer_y[2]/2) < max(st_coordinates( s)[,2]) || (maxy + buffer_y[2]/2) == max(st_coordinates( s)[, 2]))){
    int_buffer = FALSE
    message("Warning: Internal buffer has been set to FALSE because overlaps with location polygons")
  }

  if(int_buffer){
    xnodb = seq((minx - buffer_x[1]/2), (maxx + buffer_x[2]/2), by = (( (maxx + buffer_x[2]/2)- (minx - buffer_x[1]/2) )/(nnod_side_x-1)) )
    ynodb = seq((miny - buffer_y[1]/2), (maxy + buffer_y[2]/2), by = (( (maxy + buffer_y[2]/2)- (miny - buffer_y[2]/2) )/(nnod_side_y-1)) )

    edge1b = matrix(NA, ncol = 2, nrow = (nnod_side_y - 1) )
    edge1b[,1] = (minx - buffer_x[1]/2 )
    edge1b[,2] = sort(ynodb, decreasing = TRUE)[-nnod_side_y]

    edge2b = matrix(NA, ncol = 2, nrow = (nnod_side_x-1))
    edge2b[,1] = xnodb[-nnod_side_x]
    edge2b[,2] = (miny - buffer_y[1]/2)

    edge3b = matrix(NA, ncol = 2, nrow = (nnod_side_y-1))
    edge3b[,1] = (maxx + buffer_x[2]/2)
    edge3b[,2] = ynodb[-nnod_side_y]

    edge4b = matrix(NA, ncol = 2, nrow = nnod_side_x )
    edge4b[,1] = sort(xnodb, decreasing = TRUE)
    edge4b[,2] = (maxy + buffer_y[2]/2)
    polycoorb = rbind(edge1b, edge2b, edge3b, edge4b)

    border_domb = list()
    border_domb[[1]] = st_polygon(list(polycoorb))
    border_domb = st_sf(id = 2, id_location = 0, geometry = st_sfc(border_domb))

    border_dom = rbind(border_dom, border_dom_repeat, border_domb)
  }else{
    border_dom = rbind(border_dom, border_dom_repeat)
  }

  domp = rbind(border_dom, s)

  return(domp)
}
