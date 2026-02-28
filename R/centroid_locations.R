# R/centroid_locations.R

#' @noRd
#'
centroid_locations = function(domain){
  nsites = length(unique(domain[["id_location"]][which(domain[["id_location"]] > 0)]))
  location_centroid = matrix(NA, ncol = 3, nrow = nsites)

  for(i in 1:nsites){
    cpoint = sf::st_centroid(sf::st_union(domain[which(domain[["id_location"]] == i), ]))[[1]]
    location_centroid[i, ] = c(i, cpoint[1], cpoint[2])
  }

  location_centroid = data.frame(location_centroid)
  colnames(location_centroid) = c("id_location", "x_centroid", "y_centroid")

  return(location_centroid)
}
