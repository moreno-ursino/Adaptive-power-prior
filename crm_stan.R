############################################################################################
#
# Code for the usual CRM with only ESS in the PP. If ESS=0, then it is the usual CRM
#
#
############################################################################################



library(rstan)
library(ks)

# CRM with PP in stan
sm_logit <- stan_model(file="logit_regPP.stan", model_name='logit_regPP',verbose=FALSE)


# computing gamma^2 in MCMC using stan to approximate posteriors
sm_logir_gamma <- stan_model(file="logit_regN.stan", model_name='sm_logitN',verbose=FALSE)



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
  
  distancesquare=0.5*mean((sqrt(m3)-sqrt(m4))^2)*(10)
  
  return(sqrt(distancesquare))
  }
}


# necessary to compute the probabilities
invlogit <- function(x)
{
  1/(1 + exp(-x))
}


#### pklogitPP function to use in the sim file

crmPP <- function (y, doses, x, theta, D0, prob = 0.9, ESS, mod=2){ 

  # y       vector of toxicities
  # doses   pseudo doses for crm
  # x       vector of dose levels given to already accrued patients
  # theta   target prob, e.g. 0.33
  # D0      historical data, in the same form: y0, doses0, x0
  # prob    threshold of p1 to stop the trial
  # ESS     a vector with length(x) postion
  # mod     1: use the posterior mean of probabilities as estimate, 
  #         2: use the posterior mean of beta to estimate probabilities

  options = list(nchains = 4, niter = 2000, nadapt = 0.85)
  
  Npat <- length(x)
  N0 <- length(D0$x0)
  y0 = D0$y0
  doses0 = D0$doses0
  x0 = D0$x0
  
  #### ESS
  
  alphaESS <<- ESS[Npat] / N0
  
  
  ##### commensurability
  
  D = list(N=Npat, x=x, y=y, doses=doses, J=length(doses))
  


  if (Npat>100) {
    gamma1 <- commens_logit(D, D0) # only after 10 patients
  } else gamma1 = 0
  
  
  alpha0s = alphaESS*(1-gamma1)
  
  # de-comment for full power prior
  #alpha0s1 = 1

  if (Npat==1) { 
    data_s <<- list(N = Npat, y = array(y, dim=1), doses = doses, J=length(doses), x= array(x, dim=1), 
                    N0 = N0, y0=y0, doses0 = doses0, J0=length(doses0), x0=x0,
                    alpha0=alpha0s)
    
  } else {
  data_s <<- list(N = Npat, y = y, doses = doses, J=length(doses), x=x, 
                 N0 = N0, y0=y0, doses0 = doses0, J0=length(doses0), x0=x0,
                 alpha0=alpha0s)
  }
  reg <- sampling(sm_logit, data = data_s, iter = options$niter, 
                   chains = options$nchains, control = list(adapt_delta = options$nadapt))
  beta_mean = get_posterior_mean(reg, pars=c("beta"))[options$nchains+1]
  pstim_mean = get_posterior_mean(reg, pars=c("pstim"))[,options$nchains+1]
  #beta_post <- extract(reg)$beta   
  pstim_post <- extract(reg)$pstim

  
  #### computation of the probability of tox at each dose 
  if (mod == 1 ){ 
    pstim <- pstim_mean
  } else {
    pstim = invlogit(a + exp(beta_mean)*doses)
  }
  
  
  # stopping rules

  pstop <- sum(pstim_post[,1]>theta)/length(pstim_post[,1])
  stoptrial <- (pstop >= 0.9)

    
  if (stoptrial) {
    newDose = NA
  } else {
    newDose <- order((abs(pstim - theta)))[1]
  }
  parameters <- c(beta)
  names(parameters) <- c("beta")
  list(newDose = newDose, pstim = pstim, parameters = parameters,
       gamma=gamma, alphaPP = alpha0s)
}


