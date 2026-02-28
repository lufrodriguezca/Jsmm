# R/jsmm_add_triangulation.R

#' Add triangulation
#'
#' Create a domain triangulation from a polygonal object from the class sf.
#'
#' @param domain polygonal object from the class sf.
#' @param max_t_area Maximum triangle area.
#' @param min_t_angle Minimum triangle angle in degrees.
#' @return A named list with the triangulation and its corresponding polygonal sf object.
#' @export
#'
jsmm_add_triangulation = function(domain,  max_t_area, min_t_angle){
  lands_pslg = convert_to_pslg(domain)
  lands_t    = RTriangle::triangulate(p = lands_pslg, a =  max_t_area, q = min_t_angle)

  triangulation = list()
  triangulation[["ele"]]  = lands_t[["T"]]
  colnames(triangulation[["ele"]]) = c("node_1","node_2","node_3")
  triangulation[["node"]] = lands_t[["P"]]
  colnames(triangulation[["node"]]) = c("node_x", "node_y")

  nele = nrow(triangulation[["ele"]])

  dom_polygon = list()

  for(e in 1:nele){
    ele_r = triangulation[["node"]][c(unlist(triangulation[["ele"]][e, ])), ]
    dom_polygon[[e]] = sf::st_polygon(list(as.matrix(rbind(ele_r, ele_r[1, ]))))
  }

  dom_polygon = sf::st_sf(id = rep(0, nele), id_location = rep(0, nele), geometry = sf::st_sfc(dom_polygon))
  dom_polygon = sf::st_set_crs(dom_polygon, sf::st_crs(domain))

  options(warn = -1)
  c_points = sf::st_centroid(dom_polygon)
  options(warn = 0)

  p   = domain[which(domain[["id"]] != 1 ), ]

  pid = sf::st_covered_by(c_points, p)

  for(i in 1:nele){
    if(length(pid[[i]]) == 0){
      stop("Review your .shp file definition")
    }
  }

  for(i in 1:nele){
    pid_len = length(pid[[i]])
    dom_polygon[["id"]][i] = p[["id"]][pid[[i]]][[pid_len]]
  }

  p = domain[which(domain[["id"]] != 1 & domain[["id_location"]] > 0), ]
  pid = sf::st_covered_by(c_points, p)

  dom_polygon[["id_location"]] = rep(0, nele)

  for(i in 1:nele){
    v = unlist(pid[i])
    if(length(v) > 0){
      dom_polygon[["id_location"]][i] = p[["id_location"]][v]
    }
  }

  domain = list()
  domain[["polygon"]] = dom_polygon
  domain[["triangulation"]] = triangulation

  return(domain)
}


