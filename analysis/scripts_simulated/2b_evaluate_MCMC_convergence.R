library(Jsmm)
library(coda)

for(method in c("ICP", "CCP")){
  mname = paste0("simulated_", method)
  thin = 2
  samples = 250
  nChains = 5
  load(paste0(paste0("models/fitted_model_", mname, "_samples_", samples, "_thin_", thin, "_nChains_", nChains, ".RData")))
  
  mpost = convert_to_coda_object(m)
  
  #pdf(paste0("results/MCMC_convergence_", mname, ".pdf"))
  plot(mpost$Beta)
  plot(mpost$Gamma)
  if(!is.null(m$C)) plot(mpost$rho)
  summary(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf)
  summary(gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf)
  summary(gelman.diag(mpost$rho, multivariate = FALSE)$psrf)
  #dev.off()
}

bad.chains = matrix(FALSE, nrow = m$ns, ncol = 5)
for(sp in 1:m$ns){
  li = matrix(ncol = 5, nrow = 250)
  for(i in 1:5){
    for(j in 1:250){
      li[j, i] = m$postList[[i]][[j]]$li[sp]
    }
  }
  plot(NULL, xlim = c(0, 250), ylim = c(min(li), max(li)), main = sp)
  for(i in 1:5) lines(li[, i])
  me = colMeans(li)
  sd = sqrt(diag(var(li)))
  bad.chains[sp, me < max(me - 3*sd)] = TRUE
}

mean(bad.chains)
