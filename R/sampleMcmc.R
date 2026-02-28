# R/sampleMcmc.R

#' Sample MCMC
#'
#' Performs all the computations for the parameter estimation for the the Joint species movement model (Jsmm) via Markov Chain Monte Carlo (MCMC).
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @param samples Integer parameter representing the number of samples drawn from the MArkov Chain.
#' @param transient Integer parameter representing the initial part of the chain  considered as transient.
#' @param thin Integer parameter representing the thinning.
#' @param max_dt Optional positive integer value representing the number of time partitions in which the PDE should be solved.
#' @param init_pars Initial parameters for fitting the model.
#' @param nChains Number of MCMC chains to simulate.
#' @param no.data Boolean parameter indicating presence or absence of data.
#'
#' @return This functions adds to the named list m, the posterior data, thin, samples, and transient information.
#' @export
#'
sampleMcmc = function(m, samples, transient, thin, max_dt = NULL, init_pars = NULL, nChains = 1,
                      no.data = FALSE){
  update_theta = TRUE
  tpm0 = proc.time() #initial time
  if(is.null(max_dt)) max_dt = (max(m[["times"]]) - min(m[["times"]]))/50
  m = merge_rc_data(m)
  ns = m[["ns"]]
  nc = m[["nc"]]
  np = sum(nc)
  includePhylogeny = !is.null(m[["C"]])
  if(includePhylogeny) CC = m[["C"]]
  TT = m[["X"]][["trait"]]
  nt = ncol(TT)
  prior = m[["prior"]]

  muZT = prior[["muZ"]]
  VZT  = prior[["VZ"]]

  PSI = prior[["PSI"]]
  nu  = prior[["nu"]]

  if(includePhylogeny){
    rhos = prior[["rho"]][, 1]
    priorRho = prior[["rho"]][, 2]
    nr = length(rhos)
    lpriorRho = log(priorRho)

    iDDs = array(dim = c(ns,ns,nr))
    ldetDDs = numeric(nr)

    for(i in 1:nr){
      DD = rhos[i]*CC + (1 - rhos[i])*diag(ns)
      iDDs[, , i]  = solve(DD)
      ldetDDs[i] = log(det(as.matrix(DD)))
    }
  }

  XT = kronecker(TT, diag(np))

  postList = list()
  for(chain in 1:nChains){
    attempts = 0
    failures = 0
    cont = TRUE

    while(cont){
      if(attempts == 0 & !is.null(init_pars)){
        pars = init_pars
      } else {
        pars = simulate_from_prior(m)
      }
      if(no.data){
        li1 = rep(0, m[["ns"]])
      } else {
        linb = compute_likelihood(m,pars[["Beta"]], max_dt)
        li1 = linb[["log.likelihood.sp"]]
      }
      attempts = attempts + 1
      if(is.finite(sum(li1))){
        cont=FALSE
      } else {
        failures = failures + 1
      }
    }

    THETA = c()

    for(i in 1:length(nc)){
      THETA = cbind(THETA, pars[["Beta"]][[names(nc)[i]]])
    }

    SI = pars[["SI"]]
    iSI = solve(SI)
    ZT  = pars[["Gamma"]]

    if(includePhylogeny){
      RI  = round(nr/2)
      RH  = rhos[RI]
      iDD = iDDs[, , RI]
      ldetDD = ldetDDs[RI]
    } else {
      iDD = diag(ns)
      RH  = 1
    }

    ac    = array(0, c(ns, np, 2))
    kk    = array(1, c(ns, np))
    acr   = array(NA, c(ns, np))
    ns1s2 = array(0,c(ns))
    s1    = array(0,c(ns, np))
    s2    = array(0,c(ns, np, np))
    la    = array(1, c(ns, np))
    vect  = array(0,c(ns, np, np))

    for (k in 1:ns){
      vect[k, , ] = diag(np)
    }

    post = list()

    M   = TT %*% ZT
    RES = as.numeric(t(M)) - as.numeric(t(THETA))
    li2 = -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES

    print("initial likelihood: ")
    print(c(li1, li2))
    print("sampling starts")

    for (i in 1:(samples*thin + transient)){
      print(paste0("sample = ",i," out of ",samples*thin + transient,", chain = ",chain," out of ", nChains))

      if(update_theta){
        for (l in 1:np){
          cont = TRUE
          NTHETA = THETA
          for(k in 1:ns){
            nTHETA = THETA[k, ]
            mult   = stats::rnorm(1, mean = 0, sd = (kk[k, l]*sqrt(la[k, l])))
            for (l2 in 1:np){
              nTHETA[l2] = nTHETA[l2] + mult*vect[k,l2,l]
            }
            NTHETA[k, ] = nTHETA
          }

          NBeta = list()
          cnc = c(0,cumsum(nc))
          for(j in 1:length(nc)){
            if(nc[j] > 0){
              NBeta[[names(nc)[j]]] = matrix(NTHETA[, (cnc[j] + 1):cnc[j + 1]], nrow = ns)
            }
          }

          if(no.data){
            nli1 = rep(0, m[["ns"]])
          } else {
            linb = compute_likelihood(m, NBeta, max_dt)
            nli1 = linb[["log.likelihood.sp"]]
          }

          for (k in 1:ns){
            N2THETA = THETA
            N2THETA[k, ] = NTHETA[k, ]

            RES  = as.numeric(t(M)) - as.numeric(t(N2THETA))
            nli2 = -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES

            ac[k, l, 1] = ac[k, l, 1] + 1
            attempts = attempts + 1
            if(is.finite(nli1[k]) & is.finite(nli2)){
              if(runif(1) < exp((nli1[k] - li1[k] + nli2 - li2))){
                THETA[k, ] = NTHETA[k, ]
                li1[k] = nli1[k]
                li2 = nli2
                ac[k, l, 2] = ac[k, l, 2] + 1
              }
            } else {
              failures = failures + 1
            }
          }
        }
      }

      XTHETA = as.vector(t(THETA))
      iXSI   = kronecker(iDD, iSI)

      Vs  = solve(solve(VZT) + t(XT) %*% iXSI %*% XT)
      Vs  = (Vs + t(Vs))/2
      mus = Vs%*%(solve(VZT) %*% muZT + t(XT) %*% iXSI %*% XTHETA)

      ZT = matrix(MASS::mvrnorm(1,mu = mus, Sigma = Vs), ncol = np, byrow=TRUE)
      M  = TT%*%ZT
      RES = as.numeric(t(M)) - as.numeric(t(THETA))
      li2 = -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES

      RES = THETA - M
      A   = t(RES) %*% iDD %*% RES

      PSIA = PSI + A
      PSIA = (PSIA + t(PSIA))/2
      SI  = MCMCpack::riwish((nu + ns), PSIA)
      SI  = (SI + t(SI))/2
      iSI = solve(SI)
      RES = as.numeric(t(M)) - as.numeric(t(THETA))
      li2 = -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES

      if(includePhylogeny){
        RES = as.numeric(t(M)) - as.numeric(t(THETA))
        likeRho = numeric(nr)
        for(ii in 1:nr){
          likeRho[ii] = (-1/2)*(np*ldetDDs[ii] + RES%*%kronecker(iDDs[, , ii], iSI) %*% RES)
        }
        postRho = lpriorRho + likeRho

        pr  = exp(postRho)/sum(exp(postRho))
        RI  = sample(seq(1:nr), size = 1, prob = pr)
        iDD = iDDs[, , RI]

        ldetDD = ldetDDs[RI]

        RH  = rhos[RI]
        RES <- as.numeric(t(M)) - as.numeric(t(THETA))
        li2 <- -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES
      }

      if (i <= transient){
        for(k in 1:ns){
          q <- 1 + exp(-i/500)
          w <- 1 - 0.1*exp(-i/500)

          acr[k, ] = ac[k, , 2]/ac[k, , 1]
          kk[k, ] = sapply(kk[k, ] * q^(acr[k,] - 0.44), trunca)
          s1[k, ] = s1[k, ] + w*THETA[k, ]
          s2[k, , ] = s2[k, , ] + w*(THETA[k, ] %*% t(THETA[k, ]))
          ns1s2[k] = ns1s2[k] + w

          if (i > 50){
            cov = (s2[k, , ] - (s1[k, ] %*% t(s1[k, ])/ns1s2[k]))/(ns1s2[k] - 1)
            met = cov + 10^(-5)*diag(np)

            lavect      = eigen(met)
            la[k, ]     = abs(lavect$values)
            vect[k, , ] = lavect$vectors
          }
        }
        ac = ac*w
      }
      if(i > transient) {
        ii = i - transient
        if((ii %% thin) == 0){
          ii = round(ii/thin)
          pars = list()
          Beta = list()
          cnc = c(0,cumsum(nc))
          for(j in 1:length(nc)){
            if(nc[j] > 0){
              Beta[[names(nc)[j]]] = matrix(THETA[, (cnc[j] + 1):cnc[j + 1]], nrow = ns)
            }
          }
          pars[["Beta"]]  = Beta
          pars[["rho"]]   = RH
          pars[["SI"]]    = SI
          pars[["Gamma"]] = matrix(ZT, ncol = sum(nc))
          pars[["li"]]    = li1
          post[[ii]]      = pars
        }
      }
    }
    postList[[chain]] = post
  }

  tpm1 = proc.time()

  dtpm = tpm1[3] - tpm0[3]

  failure.ratio = failures/attempts

  m[["postList"]]  = postList
  m[["thin"]]      = thin
  m[["samples"]]   = samples
  m[["transient"]] = transient

  return(m)
}
