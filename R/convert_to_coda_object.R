# R/convert_to_coda_object.R

#' Convert to Coda object
#'
#' Converts the posterior data generated after fitting the jsmm model into an object to the class mcmc from the Coda R-package.
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (jsmm) analyses.
#' @param start Start of the MCMC
#' @param Beta Optional logical data type. Set to TRUE to include Beta parameters, FALSE otherwise.
#' @param Gamma Optional logical data type. Set to TRUE to include Gamma parameters, FALSE otherwise.
#' @param rho Optional logical data type. Set to TRUE to include rho parameters, FALSE otherwise.
#' @param SI Optional logical data type. Set to TRUE to include Sigma parameters, FALSE otherwise.
#' @param li Optional logical data type. Set to TRUE to include Markov chain Monte Carlo chains, FALSE otherwise.
#'
#' @return A named list with objects from the class mcmc from coda.
#' @export
#'
convert_to_coda_object = function(m, start = 1,
                               Beta = TRUE, Gamma = TRUE, rho = TRUE, SI = TRUE, li = TRUE)
{
  nChains = length(m[["postList"]])
  ns      = m[["ns"]]
  thin    = m[["thin"]]
  samples = m[["samples"]]
  start1  = m[["transient"]] + start*thin
  end1    = m[["transient"]] + samples*thin

  postListAll = m[["postList"]]

  postBeta  = list()
  postGamma = list()
  postSI    = list()
  postRho   = list()
  postLi    = list()

  covNames = m$parNames
  nc = length(covNames)
  spNames = m$sp_names
  trNames = colnames(m$X[["trait"]])
  nt = length(trNames)

  for (chain in 1:nChains){
    postList = postListAll[[chain]][start:length(postListAll[[chain]])]

    if(Beta){
      tmp = do.call(rbind, lapply(postList, function(a) as.vector(unlist(a$Beta))))
      colnames(tmp) = sprintf("B[%s, %s]", covNames[rep(1:nc, each = length(spNames))],
                              spNames[rep(1:length(spNames), nc)])
      postBeta[[chain]] = mcmc(tmp, thin = thin, start = start1, end = end1)
    }

    if (Gamma){
      tmp = do.call(rbind, lapply(postList, function(a) as.vector(a$Gamma)))
      colnames(tmp) = sprintf("G[%s, %s]", covNames[rep(1:nc, nt)], trNames[rep(1:nt, each = nc)])
      colnames(tmp) = sprintf("G[%s, %s]", trNames[rep(1:nt, nc)], covNames[rep(1:nc, each = nt)])
      postGamma[[chain]] = mcmc(tmp, thin = thin, start = start1, end = end1)
    }

    if (SI){
      tmp = do.call(rbind, lapply(postList, function(a) as.vector(a$SI) ))
      colnames(tmp) = sprintf("V[%s, %s]",covNames[rep(1:nc, nc)], covNames[rep(1:nc, each = nc)])
      postSI[[chain]] = mcmc(tmp, thin = thin, start = start1, end = end1)
    }

    if (rho){
      postRho[[chain]] = mcmc(unlist(lapply(postList, function(a) a$rho)), thin = thin,
                              start = start1, end = end1)
    }

    if (li){
      tmp = do.call(rbind, lapply(postList, function(a) a$li ))
      colnames(tmp) = sprintf("V[%s]", spNames)
      postLi[[chain]] = mcmc(tmp, thin = thin, start = start1, end = end1)
    }
  }

  mpost = list()

  if (Beta){
    mpost$Beta = mcmc.list(postBeta)
  }
  if (Gamma){
    mpost$Gamma = mcmc.list(postGamma)
  }
  if (rho){
    mpost$rho = mcmc.list(postRho)
  }
  if (SI){
    mpost$SI = mcmc.list(postSI)
  }
  if (li){
    mpost$li = mcmc.list(postLi)
  }
  return(mpost)
}
