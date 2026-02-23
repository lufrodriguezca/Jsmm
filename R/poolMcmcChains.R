# R/poolMcmcChains.R

#' Pool MCMC chains
#'
#'This function sets up the posterior data in suitable format. Used for showing posterior vs priors.
#'
#' @param postList List with the model posterior data
#' @param chainIndex Index of the MCMC chain.
#' @param start Starting point MCMC chain.
#' @param thin Thin MCMC chain
#'
#' @return A list with the formatted posterior data.
#' @export
#'

poolMcmcChains = function(postList, chainIndex = 1:length(postList),
                          start = 1, thin = 1)
{
  post = list()
  for(i in chainIndex){
    chain = postList[[i]]
    ind = seq(start, length(chain), thin)
    post = c(post, chain[ind])
  }
  return(post)
}
