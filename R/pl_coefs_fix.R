# R/pl_coefs_fix.R

#' @noRd
#'
pl_coefs_fix = function(x, el){
  pqy = x[el[, 1], 2] - x[el[, 2], 2]
  pqx = x[el[, 1], 1] - x[el[, 2], 1]
  rqx = x[el[, 3], 1] - x[el[, 2], 1]
  rqy = x[el[, 3], 2] - x[el[, 2], 2]
  c = pqx*rqy - pqy*rqx

  fix_vals = matrix(NA, ncol = 5, nrow = nrow(el))
  fix_vals[, 1] = pqx
  fix_vals[, 2] = pqy
  fix_vals[, 3] = rqx
  fix_vals[, 4] = rqy
  fix_vals[, 5] = c

  fix_vals
}
