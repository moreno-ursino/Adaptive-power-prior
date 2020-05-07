###############################################################################
###############################################################################
##
## Script for data analysis
##
## scen       scenario name (for scenario 6, see second part of the code)
##
###############################################################################
###############################################################################


vuoto <- c( " ", " ", " ", " ", " ", " ", " ", " ", " ", " ") 


scen=1

MTD_res = NULL
alloc_res = NULL
DLTs_res = NULL
#for (m in 1:30){
for (m in c(1,2,61,8,10,12,22,28,11,51)){
  #PCS
  eval(parse(text = paste("res_scen=read.table(\"MTD_mod",m,"scen",scen,".txt\")", sep="")))

  res_scen=res_scen$x
  res_scen = factor(res_scen, levels=1:6)
  MTD_res = round(rbind(MTD_res, (table(res_scen)/length(res_scen))*100))
  
  #DLT numbers
  eval(parse(text = paste("res_scen_tox=read.table(\"tox_mod",m,"scen",scen,".txt\")", sep="")))
  nDLTs = apply(as.matrix(res_scen_tox), 1, sum)
  DLTs_res = c(DLTs_res, paste(median(nDLTs), " (", quantile(nDLTs, probs = c(0.25)), ", ", quantile(nDLTs, probs = c(0.75)), ")", sep="" ))
}

library(xtable)


tabella= xtable(cbind(MTD_res, DLTs_res))
rownames( tabella ) <- c( "P_NI", "P_ESS(10)", "P_ESS(30)", "AP_L", "AP_S", "AP_MIX(0.5)", "AP_SOC1", "AP_SOC2", "AP_EB", "BCRM") 
tabella

###############
###############
###############

scen=6

MTD_res = NULL
alloc_res = NULL
DLTs_res = NULL
#for (m in 1:30){
for (m in c(1,2,61,8,10,12,22,28,11,51)){
  #PCS
  eval(parse(text = paste("res_scen=read.table(\"MTD_mod",m,"scen",scen,".txt\")", sep="")))
#res_scen=read.table("MTD_mod1scen1.txt")
res_scen=res_scen$x
if (anyNA(res_scen)) { res_scen[is.na(res_scen)] = rep(0, sum(is.na(res_scen))) }
res_scen = factor(res_scen, levels=0:6)
MTD_res = round(rbind(MTD_res, (table(res_scen)/length(res_scen))*100))


#DLT numbers
eval(parse(text = paste("res_scen_tox=read.table(\"tox_mod",m,"scen",scen,".txt\")", sep="")))
nDLTs = apply(as.matrix(res_scen_tox), 1, sum)
DLTs_res = c(DLTs_res, paste(median(nDLTs, na.rm = TRUE), " (", quantile(nDLTs, probs = c(0.25), na.rm = TRUE), ", ", 
                             quantile(nDLTs, probs = c(0.75), na.rm = TRUE), ")", sep="" ))
}

library(xtable)



tabella2= xtable(cbind(MTD_res, alloc_res, DLTs_res), digits = 2)
rownames( tabella2 ) <- c( "P_NI", "P_ESS(10)", "P_ESS(30)", "AP_L", "AP_S", "AP_MIX(0.5)", "AP_SOC1", "AP_SOC2", "AP_EB", "BCRM") 
tabella2


xtable(cbind(tabella, vuoto, tabella2))
