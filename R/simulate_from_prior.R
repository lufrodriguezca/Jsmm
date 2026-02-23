# R/simulate_from_prior.R

#' Simulate movement and observation parameters from priors.
#'
#'This function generates the movement and observation parameters from the priors rho (phylogenetic signal parameter), parameter Sigma for generating the variance-covariance matrix, and vector Gamma.
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#'
#' @return A named list with the movement and observation parameters, and the corresponding rho, Sigma, and Gamma priors.
#' @export
#'
simulate_from_prior = function(m){

  prior = m[["prior"]]
  muZ = prior[["muZ"]]
  VZ  = prior[["VZ"]]

  PSI = prior[["PSI"]]
  nu  = prior[["nu"]]

  includePhylogeny = !is.null(m[["C"]])
  if(includePhylogeny){
    rhos = prior[["rho"]][, 1]
    prior_rho    = prior[["rho"]][, 2]
  }

  np = sum(m[["nc"]])
  Ins = diag(m[["ns"]])

  if(!is.null(m[["C"]])){
    rho = sample(size = 1, prob = prior_rho, x = rhos)
  }else{
    rho = 0
  }

  if(is.null(m[["C"]])){
    C = diag(m[["ns"]])
  }else{
    C = m[["C"]]
  }

  SI = riwish(v = nu, S = PSI)

  tr = m[["X"]][["trait"]]

  nt = dim(tr)[2]

  Z = matrix(mvrnorm(1, mu = muZ, Sigma = VZ), ncol = np, byrow = TRUE)

  M = tr%*%Z
  em = matrix(t(M), nrow = 1)

  THETAS = matrix(mvrnorm(1, mu = em, Sigma = kronecker(rho*C + (1 - rho)*Ins, SI)), nrow = 1)
  THETAS = matrix(THETAS, ncol = np, byrow = TRUE)

  cs = cumsum(m$nc)
  vars = names(cs)
  cs = c(0, cs)

  Beta = list()

  for(i in seq_len(length(vars))){
    if(cs[i + 1] > cs[i]) Beta[[vars[[i]]]] = matrix(THETAS[, (cs[i] + 1):cs[i + 1]], nrow = m[["ns"]])
  }

  pars = list()

  pars[["Beta"]]   = Beta
  pars[["rho"]]    = rho
  pars[["SI"]]     = SI
  pars[["Gamma"]]  = Z

  return(pars)
}


