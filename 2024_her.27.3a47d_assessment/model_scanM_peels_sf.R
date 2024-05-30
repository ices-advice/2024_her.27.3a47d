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

runName     <- 'NSAS_HAWG2024_sf_scanM_peels'

SMSkeyRun <- 2023

source('utilities_model.R')
source('utilities_model_config.R')

nPeels <- 11

### ============================================================================
### Construct packages object
### ============================================================================

packages <- as.data.frame(taf.session()$packages)

save(packages,
     file=file.path(model.save,paste0("packages_",runName,".RData")))

### ============================================================================
### loop on peels
### ============================================================================

NSH.sams  <- new("FLSAMs")
NSHs      <- new("FLStocks")
NSH.ctrls <- list()

for(iPeel in 1:nPeels){
  
  ### ============================================================================
  ### load object
  ### ============================================================================
  load(file.path(data.save,paste0("data_sms2023_sf.RData")))  # NSH, NSH.tun
  
  endYear <- range(NSH)['maxyear']
  runNameLoop <- ac(endYear-(iPeel-1))
  
  NSH     <- window(NSH, end=endYear-(iPeel-1))
  NSH.tun$HERAS <- window(NSH.tun$HERAS, end=endYear-(iPeel-1))
  NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=endYear-(iPeel-1))
  NSH.tun$`LAI-ORSH` <- window(NSH.tun$`LAI-ORSH`, end=endYear-(iPeel-1))
  NSH.tun$`LAI-BUN` <- window(NSH.tun$`LAI-BUN`, end=endYear-(iPeel-1))
  NSH.tun$`LAI-CNS` <- window(NSH.tun$`LAI-CNS`, end=endYear-(iPeel-1))
  NSH.tun$`LAI-SNS` <- window(NSH.tun$`LAI-SNS`, end=endYear-(iPeel-1))
  NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=endYear-(iPeel-1))
  NSH.tun$`IBTS-Q1` <- window(NSH.tun$`IBTS-Q1`, end=endYear-(iPeel-1)+1)
  NSH.tun$IBTS0     <- window(NSH.tun$IBTS0, end=endYear-(iPeel-1)+1)
  
  pg <- range(NSH)[2]

  ### ============================================================================
  ### Construct control object
  ### ============================================================================
  
  NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)

  NSH.ctrls[[ac(runNameLoop)]] <- NSH.ctrl
  
  ### ============================================================================
  ### Run model
  ### ============================================================================
  
  # pre-run assessment for initial values
  
  mOrig <- NSH@m
  inc <- 0.01
  startSeq  <- -0.3
  endSeq    <- 0.3
  addM <- seq(startSeq,endSeq,inc)
  
  res <- run_scanM(addM,NSH,NSH.tun,NSH.ctrl,mOrig)
  
  save(NSH, NSH.tun, NSH.ctrl,res,
       file=file.path(model.save,
                      paste0(runName,'_peel',endYear-(iPeel-1),"_sf.RData")))
  
  # store optimal result
  df.nlogl       <- as.data.frame(unlist(lapply(res$NSH.sams,nlogl)))
  colnames(df.nlogl) <- 'nlogl'
  df.nlogl$addM  <- rownames(df.nlogl)
  
  idxOpt <- which(min(df.nlogl$nlogl) == df.nlogl$nlogl)
  
  NSH.sam.opt   <- res$NSH.sams[[ac(df.nlogl$addM[idxOpt])]]
  M             <- res$NSH.m[[ac(df.nlogl$addM[idxOpt])]]
  
  NSH.opt           <- NSH
  NSH.opt@stock.n   <- NSH.sam.opt@stock.n[,ac(range(NSH.opt)["minyear"]:range(NSH.opt)["maxyear"])]
  NSH.opt@harvest   <- NSH.sam.opt@harvest[,ac(range(NSH.opt)["minyear"]:range(NSH.opt)["maxyear"])]
  NSH.opt@m         <- M
  
  NSH.sams[[ac(runNameLoop)]]  <- NSH.sam.opt
  NSHs[[ac(runNameLoop)]]      <- NSH.opt
}

save(NSH.ctrls,
     file=file.path(model.save,
                    paste0("config_",runName,"_sf.RData")))

save(NSHs, NSH.tun, NSH.ctrls, NSH.sams,
     file=file.path(model.save,
                    paste0(runName,".RData")))