# R/compute_integral.R

#' @noRd
#'
compute_integral = function(m, XBeta){
  integral = assemble_integral(k_v = c(XBeta[["habitat_preference"]]) , detB = m[["f_data"]][["det_B"]],
                    nnod = nrow(m[["domain"]][["triangulation"]][["node"]]), ns = m[["ns"]],
                    nele = nrow(m[["domain"]][["triangulation"]][["ele"]]), el_v = c(m[["domain"]][["triangulation"]][["ele"]]))

  integral = matrix(integral, ncol = m[["ns"]])

  return(integral)
}
