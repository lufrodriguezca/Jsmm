# R/compute_conversion_factors.R

#' Conversion factors computation
#'
#' @description This function implements the conversion factor \eqn{Z^* = \frac{1}{Ap}} method described in Rodriguez, Gårdman, et al. (2026).  \eqn{Z^*} represents the population density that results in the expected number of trapped individuals being \eqn{E[N] = 1}. \eqn{A} and \eqn{p} represent the domain area and probability by which the individual would become trapped according to the continuous capture process, respectively.
#'
#' @references
#' Rodriguez L.F., Gårdman V., Roslin T., & Ovaskainen O. (2026). Malaise trap samples of 1000 individuals per week suggest 4 million insects per hectare in the boreal zone.
#'
#' @param m Jsmm named list that contains the fitted model.
#' @param time1 positive real number that defines the starting time for the time duration period T.
#' @param time2 positive real number that defines the end time for the time duration period T.
#' @param max_dt  positive real number that describes the maximal time resolution used in numerical computations.
#'
#' @return A matrix with the computations for the species defined in m.
#' @export
#'
compute_conversion_factors = function(m, max_dt = NULL, time1 = 0, time2 = 1){
  max_dt = 0.0205
  delta  = time2 - time1
  ltimes = seq(from = time1, to = time2, length.out = 1 + ceiling(delta/max_dt))
  post = m[["postList"]]

  number.per.hectare = c()

  for(chain in 1:length(post)){
    for(sa in 1:length(post[[chain]])){
      print(c(chain,sa,length(post[[chain]])))
      pars = post[[chain]][[sa]]
      pars[["Beta"]][["mortality"]] = -100 + 0*pars[["Beta"]][["mortality"]]
      pre = precomputation_KM(m)
      nnod = nrow(m[["domain"]][["triangulation"]][["node"]])

      v = matrix(0, nrow = nnod, ncol = m[["ns"]])

      XBeta_start = multiply_XBeta(m = m, Beta = pars[["Beta"]], i = time1)
      XBeta_end   = multiply_XBeta(m = m, Beta = pars[["Beta"]], i = time2)

      #RELEASE INDIVIDUALS AT UNIFORM DENSITY
      XBeta = XBeta_start
      integral = compute_integral(m = m, XBeta = XBeta)
      for(sp in 1:m[["ns"]]){
        v[, sp] = 1/sum(integral[, sp])
      }

      sum(v[, 1]*integral[, 1])

      #INNER LOOPS STARTS
      for(i2 in 1:(length(ltimes) - 1)){
        #COMPUTE D1 AND D2
        current_time = (ltimes[i2] + ltimes[i2 + 1])/2
        re.dt = ltimes[i2 + 1] - ltimes[i2]
        w2 = (current_time - time1)/(time2 - time1)
        w1 = 1 - w2
        XBeta = lapply(XBeta_start,"*", w1)
        tmp = lapply(XBeta_end, "*", w2)

        for(i in 1:length(XBeta)) XBeta[[i]] = XBeta[[i]] + tmp[[i]]
        K = compute_Stiffness(m = m, XBeta = XBeta, pre)
        M = compute_Mass(m = m, XBeta = XBeta, pre)
        integral = compute_integral(m = m, XBeta = XBeta)
        D1 = list()
        D2 = list()

        for(sp in 1:m[["ns"]]){
          D1[[sp]] = t(M[[sp]]) - (2/3)*re.dt*t(K[[sp]])
          D2[[sp]] = t(M[[sp]]) + (1-(2/3))*re.dt*t(K[[sp]])
          u = D2[[sp]]%*%v[,sp]
          v[, sp] = as.numeric(Matrix::solve(D1[[sp]], u))
        }
      }

      integral = compute_integral(m = m, XBeta = XBeta)
      sum(v[, 1]*integral[, 1])
      pr.not.trapped = rep(NA, m$ns)
      for(i in 1:m$ns){pr.not.trapped[i] = sum(v[, i]*integral[, i])}
      pr.trapped = 1 - pr.not.trapped
      #probability of the individual becoming trapped

      e.per.hectare = pr.trapped * sum(st_area(m[["domain"]][["polygon"]][["geometry"]]))/(100^2)
      #expected number of trapped individuals, assuming one individual per hectare

      number.per.hectare = rbind(number.per.hectare,1/e.per.hectare)
      #the number of individuals per hectare that leads to expected number of individuals captured being one
    }
  }
  colnames(number.per.hectare) = m[["sp_names"]]

  return(number.per.hectare)
}
