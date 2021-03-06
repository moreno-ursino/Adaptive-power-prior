###############################################################################
###############################################################################
##
## Script to run model estimations
##
###############################################################################
###############################################################################

  # add ESS for each person of the trial ESS[1] = ESS I want if I have only one patient in the new trial


  CI=FALSE
  cohort=1
  Ntot=N=30
  theta=target
  TR=1000
  

  seed = 15478 # or what you prefer
  

  simparallel <- function(tr, simulatedData, D0, cohort, Ntot, a, doses, theta, ESS, seed){
    
    # tr                    number of trial, usually index for the "for"
    # simulatedData         simulated TR data for actual trial
    # D0                    simulated 1 dataset of historical data
    # cohort                cohort size
    # Ntot                  maximum number of patients
    # a                     logistic intercept
    # doses                 pseudo doses for CRM
    # theta                 toxicity threshold
    # ESS                   ESS function
    # seed                  set.seed
    
    
    J <- length(doses)
    tox <- simulatedData[[tr]]           #for detail, see how we simulated these matrices

    
    
    
    
    x <- rep(1, cohort)
    y <- tox[cbind(1:length(x), x)]
    M = Ntot/cohort
    
    alphaPPtot = rep(0, cohort)     # to check all the alpha during the trial


    stage1 = FALSE       # here STAGE 1 is before the first toxicity
    for (i in 2:M) {
      j = (cohort * (i - 1) + 1):(cohort * i)
      if (stage1) {
        x <- c(x, rep(min((max(x) + 1), length(doses)), 
                      cohort))
        y <- c(y, tox[cbind(j, x[j])])

        alphaPPtot <- c(alphaPPtot, rep(0,cohort))

        if (any(y == "1")) {
          stage1 <- FALSE
        }
      } else {    
        
        set.seed(seed + tr + length(x)) 
        results <- crmPP(y=y, doses=doses, x=x, theta=theta, D0=D0, ESS=ESS) 
        if (is.na(results$newDose) == "TRUE") break
        newdose <- min(results$newDose, max(x) + 1)
        x <- c(x, rep(newdose, cohort))
        y <- c(y, tox[cbind(j, x[j])])
        alphaPPtot <- c(alphaPPtot, rep(results$alphaPP, cohort))

      }
    }
    trial <- paste("trial:", tr, sep = "")
    if (length(x) < Ntot) {
      nstop <- N - length(x)
      MTD = 0
      doseLevels =  c(x, rep(NA, nstop))
      toxicity = c(y, rep(NA, nstop))
      alphaPPtot =  c(alphaPPtot, rep(NA, nstop))
    } else {
      results = crmPP(y=y, doses=doses, x=x, theta=theta, D0=D0, ESS=ESS)
      MTD = results$newDose
      doseLevels = x
      toxicity = y
      alphaPPtot <- c(alphaPPtot, results$alphaPP)
    }
    return(list(MTD=MTD, doseLevels=doseLevels, toxicity=toxicity, alphas = alphaPPtot))
  }
  
  
  start_time1 <- Sys.time()
  
  resultats <- mclapply(1:TR, simparallel,simulatedData=simulatedData, D0=D0, 
                       cohort=cohort, Ntot=N, a=a, doses=doses, theta=theta, ESS=ESS,
                       seed= seed, mc.cores=8)
  
  end_time1 <- Sys.time()
  start_time1 - end_time1
  

  MTD_sim <- NULL    # median
  doselevels_sim <- NULL   # mean
  tox_sim <- NULL
  alpha_sim <- NULL
  
  for (tr in 1:TR){
    MTD_sim <- c(MTD_sim, resultats[[tr]]$MTD)
    doselevels_sim <- rbind(doselevels_sim, resultats[[tr]]$doseLevels)
    tox_sim <- rbind(tox_sim, resultats[[tr]]$toxicity)
    alpha_sim <- rbind(alpha_sim, resultats[[tr]]$alphas)
  }
  
  write.table(MTD_sim, file=paste("MTD_mod",modello,"scen",scen,".txt", sep="") )
  write.table(doselevels_sim, file=paste("doselevels_mod",modello,"scen",scen,".txt", sep="") )
  write.table(tox_sim, file=paste("tox_mod",modello,"scen",scen,".txt", sep="") )
  write.table(alpha_sim, file=paste("alpha_mod",modello,"scen",scen,".txt", sep="") )
  
