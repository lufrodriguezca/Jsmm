# R/set_prior_parameters.R

#' @noRd
#'
set_prior_parameters = function(m, rho = NULL, PSI = NULL, nu = NULL, muZ = NULL, VZ = NULL){
  np = sum(m[["nc"]])
  prior = list()

  if(!is.null(m[["C"]])){
    if(is.null(rho)){
      rhos = seq(0, 1, 0.01)
      nr   = length(rhos)
      prior_rho = rep(0.5/(nr - 1), nr)
      prior_rho[1] = 0.5
      rho = cbind(rhos,prior_rho)
      colnames(rho) = c("rho","probability")
    }
    prior[["rho"]] = rho
  }

  if(is.null(PSI)){
    PSI = diag(np)
  }
  prior[["PSI"]] = PSI
  if(is.null(nu)){
    nu = np + 5
  }
  prior[["nu"]] = nu

  tr = m[["X"]][["trait"]]
  nt = ncol(tr)
  if(is.null(muZ)){
    muZ = matrix(0, nrow = nt*np)
  }
  prior[["muZ"]] = muZ

  if(is.null(VZ)){
    VZ  = diag(np*nt)
  }
  prior[["VZ"]] = VZ

  m[["prior"]] = prior

  return(m)
}

