# R/compute_Mass.R

#' @noRd
#'
compute_Mass = function(m, XBeta, pre){
  detB = m[["f_data"]][["det_B"]]
  nele = nrow(m[["domain"]][["triangulation"]][["ele"]])
  MT = pre[["MT"]]
  xPos = pre[["xPos"]]
  ML = pre[["ML"]]
  vals = pre[["vals"]]

  Ms = list()

  for(sp in 1:m[["ns"]]){
    vals = vals*0
    k = XBeta[["habitat_preference"]][, sp]
    vals = assemble_mass(k = k, ML = ML, detB = detB, vals = vals, nele = nele, xPos = xPos)
    MT@x = vals
    Ms[[sp]] = MT
  }
  return(Ms)
}

