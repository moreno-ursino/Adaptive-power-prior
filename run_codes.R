###############################################################################
###############################################################################
##
## Script to run simulation. It will run all models sequentially. 
## The user can also select the piece of the script corresponding to the 
## desired method
##
###############################################################################
###############################################################################

# script that permits to generate scenarios datasets
source(file="generation_scen.R")

require(parallel)

source(file="crm_stan.R")

#P_NI
modello=1 
ESS = rep(0, 40)
source(file="prototyp_nsim.R")

#P_ESS(10)
modello=2
ESS = rep(10, 40)
source(file="prototyp_nsim.R")

#P_ESS(30)
modello=61
ESS = rep(30, 40)
source(file="prototyp_nsim.R")

#AP_L
modello=8
ESS = 1:40
source(file="crm_stan_alpha.R")
source(file="prototyp_nsim.R")

#AP_S
modello=10
ESS = 1:40
source(file="crm_stan_alpha5.R")
source(file="prototyp_nsim.R")

#AP_MIX(0.5)
modello=12
ESS = 1:40
source(file="crm_stan_rev.R")
source(file="crm_stan_mix.R")
w_mix=0.5
source(file="prototyp_nsim.R")


#AP_SOC1
modello=22
ESS = 1:40
soglia_occam=0.2
source(file="crm_stan.R")
source(file="crm_stan_alpha6.R")
source(file="crm_stan_alpha5.R")
source(file="prototyp_nsim.R")


#AP_SOC2
modello=28
ESS = pmin(rep(20,40), 1:40)
soglia_occam=0.2
source(file="prototyp_nsim.R")

#AP_EB
modello=11
ESS = 1:40  # not used but to not modify the fucntion
source(file="EB.R")
source(file="prototyp_nsim.R")

#BCRM
modello=51
ESS = pmin(rep(20,200), 1:200)
source(file="Lee_method2.R")
source(file="prototyp_nsim3.R")
