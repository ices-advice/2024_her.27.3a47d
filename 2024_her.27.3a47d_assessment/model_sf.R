## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

rm(list=ls())

library(icesTAF)

taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)

# taf.session()

library(methods)


mkdir("model")

data.save               <- file.path(".", "data")
model.save              <- file.path(".",'model')

SMSkeyRuns  <- 2023
retroFlag   <- TRUE

runName <- 'NSAS_HAWG2024_sf'

source('utilities_model_config.R')

### ============================================================================
### Construct packages object
### ============================================================================

packages <- as.data.frame(taf.session()$packages)

save(packages, 
     file=file.path(model.save,paste0("packages_",runName,".RData")))

### ============================================================================
### Construct single fleet control object
### ============================================================================

addM    <- 0.02

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun

NSH@m   <- NSH@m + addM

pg <- range(NSH)[2]

NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)
NSH.ctrl@residuals <- TRUE

### ============================================================================
### Run single fleet model
### ============================================================================

# Run model
NSH.sam     <- FLSAM(NSH, NSH.tun, NSH.ctrl)
# NSH.sim         <- simulate(NSH,NSH.tun,NSH.ctrl,n=10)
# temp <- monteCarloStock(NSH,NSH.tun,NSH.sam,10)
# 
# NSH.sim         <- simulate(NSH,NSH.tun,NSH.ctrl,n=nits)

NSH@stock.n <- NSH.sam@stock.n[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]
NSH@harvest <- NSH.sam@harvest[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]

save(NSH, NSH.tun, NSH.ctrl, NSH.sam, 
     file=file.path(model.save,paste0(runName,'.RData')))

### ============================================================================
### run single fleet retro
### ============================================================================

if(retroFlag){
  # retro
  n.retro.years <- 10  # Number of years for which to run the retrospective
  NSH.ctrl@residuals <- FALSE
  NSH.retro <- retro(NSH, NSH.tun, NSH.ctrl, retro=n.retro.years)
  
  save(NSH.retro,
       file=file.path(model.save,paste0(runName,'_retro.RData')))
}
