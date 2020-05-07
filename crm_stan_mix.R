###############################################################################
###############################################################################
##
## For mixture priors. Required: crm_stan.R 
##
###############################################################################
###############################################################################




library(rstan)


# CRM with PP in stan
sm_logit_mix <- stan_model(file="logit_regPP_mix.stan", model_name='logit_regPP',verbose=FALSE)


denominator_int <- function(theta, x, y, doses, alpha, a){
  loglikelihood <- NULL
  
  p_d=NULL
  for (j in 1:length(doses)) {
    p_d[j] = 1/(1+exp(- (a + exp(theta)*doses[j]))) 
    if (p_d[j]==0) p_d[j] = .Machine$double.xmin
    if (p_d[j]==1) p_d[j] = 0.9999998
  }  
  #print(p_d)
  
  
  for (i in 1:length(x)){
    loglikelihood[i] = y[i]*log(p_d[x[i]]) + (1-y[i])*log(1-p_d[x[i]]) 
  }
  
  lik = alpha*sum(loglikelihood)+dnorm(theta, sd=sqrt(1.34), log=TRUE)
  #lik = dnorm(theta, sd=1.34, log=TRUE)
  exp(lik)
}




#### pklogitPP function to use in the sim file

crmPP <- function (y, doses, x, theta, D0, prob = 0.9, ESS, mod=2, a=3){

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
  


  if (Npat>10) {
    gamma1 = commens_logit(D, D0) # only after 10 patients
    alpha0s = alphaESS*max((1-gamma1),0)
  } else alpha0s = 0
  
 
  if (alpha0s ==0) {
    dem=1
  } else  dem= integrate(Vectorize(denominator_int, vectorize.args = "theta"), lower=-Inf, upper = Inf, x=x0, y=y0, 
                 doses=doses0, alpha=alpha0s, a=a)$value
  
  # de-comment for full power prior
  #alpha0s1 = 1

  if (Npat==1) {
    data_s <<- list(N = Npat, y = array(y, dim=1), doses = doses, J=length(doses), x= array(x, dim=1),
    N0 = N0, y0=y0, doses0 = doses0, J0=length(doses0), x0=x0,
    alpha0=alpha0s, int_dem=dem, w_mix=w_mix)
    
  } else {
    data_s <<- list(N = Npat, y = y, doses = doses, J=length(doses), x=x,
    N0 = N0, y0=y0, doses0 = doses0, J0=length(doses0), x0=x0,
    alpha0=alpha0s, int_dem=dem, w_mix=w_mix)
  }

  reg <- sampling(sm_logit_mix, data = data_s, iter = options$niter,
                   chains = options$nchains, control = list(adapt_delta = options$nadapt))
  beta_mean = get_posterior_mean(reg, pars=c("beta"))[options$nchains+1]
  pstim_mean = get_posterior_mean(reg, pars=c("pstim"))[,options$nchains+1]
  #beta_post <- extract(reg)$beta   
  #pstim_post <- extract(reg)$pstim

  
  #### computation of the probability of tox at each dose 
  if (mod == 1 ){ 
    pstim <- pstim_mean
  } else {
    pstim = invlogit(a + exp(beta_mean)*doses)
  }
  
  
  # # stopping rules
  # pstim_sum <- matrix(0, ncol = options$nchains * options$niter/2, nrow=1)
  #                     #nrow = length(doses))
  # p_sum <- NULL
  # parmt1 = sampl1$b[, 1] + sampl1$b[, 2] * log(doses[1])
  # parmt2 = sampl1$sigma
  # for (i in 1:ncol(pstim_sum)) {
  #   pstim_sum[1, i] <- integrate(f2_logit, -Inf, Inf, lambda1 = sampl2$beta2[i], 
  #                                lambda2 = sampl2$beta3[i], parmt1 = parmt1[i], parmt2 = parmt2[i])$value
  # }
  # pstop <- checking1(pstim_sum[1, ], target = theta, error = 0)
  # stoptox <- (pstop >= prob)
  # stoptrial <- stoptox


  ## add here for CI for other doses and not only the first.
  stoptrial = FALSE
    
  if (stoptrial) {
    newDose = NA
    message("The trial stopped based on the stopping rule \n \n")
  } else {
    newDose <- order((abs(pstim - theta)))[1]
  }
  parameters <- c(beta)
  names(parameters) <- c("beta")
  list(newDose = newDose, pstim = pstim, parameters = parameters,
       gamma=gamma, alphaPP = alpha0s)
}


