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

source('utilities_model_config.R')

### ============================================================================
### ============================================================================
### sf run without 2024 IBTS-Q1 age 1
### ============================================================================
### ============================================================================

runName <- 'NSAS_HAWG2024_sf_no2024IBTSQ1'
SMSkeyRuns  <- 2023

### ============================================================================
### Construct single fleet control object
### ============================================================================

addM    <- 0.02

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun

NSH.tun$`IBTS-Q1` <- window(NSH.tun$`IBTS-Q1`, end=2023)

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
### ============================================================================
### sf run with 2019 SMS keyrun
### ============================================================================
### ============================================================================

runName <- 'NSAS_HAWG2024_sf_SMS2019'
SMSkeyRuns  <- 2019

### ============================================================================
### Construct single fleet control object
### ============================================================================

addM    <- 0.05

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
### ============================================================================
### sf retro run with 2019 SMS keyrun and addM=0.06 as for HAWG2023
### ============================================================================
### ============================================================================

runName <- 'NSAS_HAWG2024_sf_SMS2019_2'
SMSkeyRuns  <- 2019


addM    <- 0.06

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun

NSH@m   <- NSH@m + addM

pg <- range(NSH)[2]

NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)
NSH.ctrl@residuals <- TRUE

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

# retro
n.retro.years <- 10  # Number of years for which to run the retrospective
NSH.ctrl@residuals <- FALSE
NSH.retro <- retro(NSH, NSH.tun, NSH.ctrl, retro=n.retro.years)

save(NSH.retro,
     file=file.path(model.save,paste0(runName,'_retro.RData')))

### ============================================================================
### ============================================================================
### sf HAWG2023
### ============================================================================
### ============================================================================

runName <- 'NSAS_HAWG2024_sf_IBP2021'
SMSkeyRuns  <- 2019


addM    <- 0.06

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun

NSH     <- window(NSH, end=2020)
NSH.tun$HERAS <- window(NSH.tun$HERAS, end=2020)
NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=2020)
NSH.tun$`LAI-ORSH` <- window(NSH.tun$`LAI-ORSH`, end=2020)
NSH.tun$`LAI-BUN` <- window(NSH.tun$`LAI-BUN`, end=2020)
NSH.tun$`LAI-CNS` <- window(NSH.tun$`LAI-CNS`, end=2020)
NSH.tun$`LAI-SNS` <- window(NSH.tun$`LAI-SNS`, end=2020)
NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=2020)
NSH.tun$`IBTS-Q1` <- window(NSH.tun$`IBTS-Q1`, end=2020+1)
NSH.tun$IBTS0     <- window(NSH.tun$IBTS0, end=2020+1)

NSH@m   <- NSH@m + addM

pg <- range(NSH)[2]

NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)
NSH.ctrl@residuals <- TRUE

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
