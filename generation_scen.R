###############################################################################
###############################################################################
##
## Script for data simulation
##
## scen       scenario name
## TR         number of trial generated
## N          maximum sample size for each trial
## preal      generating probabilities of toxicity (vector of the same length of the dose panel)
##
###############################################################################
###############################################################################

scen=1

TR = 1000
N = 30
preal = c(0.001, 0.01, 0.05 ,0.07, 0.2, 0.4)
SED = 158


###############   sim.data function adapted

sim.data.custom<-function (preal, N, TR, seed) 
{
  
  resScen <- list()
  for (tr in 1:TR) {
    set.seed(seed + tr)
    tox <- NULL
    
    for (i in 1:N) {
      tox <- rbind(tox, rbinom(length(preal), size=1, prob=preal))
    }

    rownames(tox) <- paste("pid", 1:N)
    colnames(tox) <- paste("dose", 1:length(preal))
    resScen[[tr]] <- tox
  }
  return(resScen)
}


########################################


simulatedData <- sim.data.custom(preal, N, TR, seed=SED) 
save.image()



########################  generating historical data

library(dfcrm)

SEED = 841
target=0.2
prior=c(0.05, 0.07, 0.2, 0.4, 0.5, 0.55)

simulatedData0 <- crmsim(PI=c( 0.05, 0.07, 0.2, 0.4, 0.5, 0.55), prior=prior, target, N, mcohort=1, model="logistic", x0=1, seed=SEED)
#simulatedData0

a=3
doses <- (log(prior/(1-prior)) - a)/1

D0 = list(x0 = simulatedData0$level, y0=simulatedData0$tox, doses0= doses)

save.image()
