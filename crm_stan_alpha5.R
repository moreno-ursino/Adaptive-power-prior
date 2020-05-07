###############################################################################
###############################################################################
##
## Distance square-root. 
## Required: crm_stan.R and crm_stan_alpha.R 
##
###############################################################################
###############################################################################



#### check commensurability
commens_logit <- function(D, D0){
  
  it=2000
  
  D$NN0 = 1
  D$a = 3
  data_m2 <- D # current data, full power and non info prior (flat)  
  reg_m2 <- sampling(sm_logir_gamma, data = data_m2, iter = it, 
                     chains = 4, control = list(adapt_delta = 0.85))
  m2_sampl <- extract(reg_m2, pars=c("beta"))$beta
  if (any(summary(reg_m2)$summary[,"Rhat"] > 1.1))   reg_m2 <- sampling(sm_logir_gamma, data = data_m2, iter = it, chains = 4, control = list(adapt_delta = 0.85))
  if (any(summary(reg_m2)$summary[,"Rhat"] > 1.1))   reg_m2 <- sampling(sm_logir_gamma, data = data_m2, iter = it, chains = 4, control = list(adapt_delta = 0.85))
  
  
  data_m3 <- list(N=length(D0$x0), x=D0$x0, y=D0$y0, doses=D0$doses0, J=length(D0$doses0),
                  NN0 = D$N/length(D0$x0), a=3)
  reg_m3 <<- sampling(sm_logir_gamma, data = data_m3, iter = it, 
                     chains = 4, control = list(adapt_delta = 0.85))
  if (any(summary(reg_m3)$summary[,"Rhat"] > 1.1))   reg_m3 <- sampling(sm_logir_gamma, data = data_m3, iter = it, chains = 4, control = list(adapt_delta = 0.85))
  if (any(summary(reg_m3)$summary[,"Rhat"] > 1.1))   reg_m3 <- sampling(sm_logir_gamma, data = data_m3, iter = it, chains = 4, control = list(adapt_delta = 0.85))
  
  m3_sampl <- extract(reg_m3, pars=c("beta"))$beta
  
  if (any(summary(reg_m3)$summary[,"Rhat"] > 1.1) | any(summary(reg_m2)$summary[,"Rhat"] > 1.1)){
    return(-1)
  } else {

  npoints = 10000
  # simulate number in hypercube
  
  points_integral <- runif(npoints, min=-5, max=5)
  
  
  m3 <<- pmax(kde(x=as.numeric(m2_sampl), eval.points=points_integral)$estimate , rep(0, npoints))
  m4 <<- pmax(kde(x=as.numeric(m3_sampl), eval.points=points_integral)$estimate , rep(0, npoints))
  
  #mean(m3$estimate)*(20*10*1) 
  #mean(m4$estimate)*(20*10*1)
  
  distancesquare=0.5*mean((sqrt(m3)-sqrt(m4))^2)*(10)
  
  return(distancesquare^(1/4))
  }
}

