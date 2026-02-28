# R/affine_mapping.R

#' @noRd
#'
affine_mapping = function(triangulation){
  el = matrix(unlist(triangulation[["ele"]]), ncol = 3)
  x  = matrix(unlist(triangulation[["node"]]), ncol = 2)

  nele = nrow(el)
  nnod = nrow(x)

  r = 10^(-20)

  xr_1 = r*cos(0)
  yr_1 = r*sin(0)

  xr_2 = r*cos(pi/2)
  yr_2 = r*sin(pi/2)

  B_j    = matrix(NA, ncol = 4, nrow = nele)
  det_B  = rep(NA, nele)
  metric = matrix(NA, ncol = 4, nrow = nele)

  for(i in 1:nele){
      B_j[i,]    = c(x[el[i, 2], 1] - x[el[i, 1], 1], x[el[i, 3], 1] - x[el[i, 1], 1], x[el[i, 2], 2] - x[el[i, 1], 2], x[el[i, 3], 2] - x[el[i, 1], 2])
      det_B[i]   = B_j[i, 1]*B_j[i, 4] - B_j[i, 2]*B_j[i, 3]
      metric[i,] = c((1/det_B[i])*c(x[el[i, 3], 2] -x[el[i, 1], 2], -(x[el[i, 2], 2] -x[el[i, 1], 2]), -(x[el[i, 3], 1] -x[el[i, 1], 1] ), x[el[i, 2], 1] - x[el[i, 1], 1]))
  }

  fix_vals = pl_coefs_fix(x, el)

  pe_1 = matrix(c(c(1, 2, 3), c(2, 3, 1), c(3, 1, 2)), ncol = 3, byrow = T)

  k = 0

  for(i in 1:nele){
    for(j in 1:3){
      k = k + 1
      el_perm = el[i, pe_1[j, ]]

      B_jperm = c(x[el_perm[2], 1] - x[el_perm[1], 1],
                 x[el_perm[3], 1] - x[el_perm[1], 1],
                 x[el_perm[2], 2] - x[el_perm[1], 2],
                 x[el_perm[3], 2] - x[el_perm[1], 2])

      p1 = matrix(B_jperm, ncol = 2, byrow = T)%*%c(xr_1, yr_1) +  c(x[el_perm[1], 1], x[el_perm[1], 2])
      p2 = matrix(B_jperm, ncol = 2, byrow = T)%*%c(xr_2, yr_2) +  c(x[el_perm[1], 1], x[el_perm[1], 2])
    }
  }

  f_data = list()
  f_data[["B_j"]] = B_j
  f_data[["det_B"]] = det_B
  f_data[["metric"]] = metric
  f_data[["fix_vals"]] = fix_vals

  return(f_data)
}

