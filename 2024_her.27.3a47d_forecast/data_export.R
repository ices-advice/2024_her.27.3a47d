### ============================================================================
### Export data
### ============================================================================

#rm(list=ls())

setwd("data")

write.taf(as.data.frame(unlist(referencePoints)), "refpoints.csv")
write.taf(catchAreaTab, "NSAS_WBSS_catch_split_uptake.csv")
write.taf(transferTab, "C_fleet_transfer.csv")
write.taf(TACTab, "TAC_effective.csv")
write.taf(adviceTab, "advise_TAC.csv")

setwd("..")
