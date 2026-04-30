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

runName     <- 'NSAS_HAWG2024_sf_scanM_smsAll'

#mkdir(file.path("model","assessment",runName))
#model.assessment.save.run   <- file.path(".",'model',"assessment",runName)

SMSkeyRunsMat <- c(2010,2013,2016,2019,2023)

source('utilities_model.R')
source('utilities_model_config.R')

### ============================================================================
### Construct packages object
### ============================================================================

packages <- as.data.frame(taf.session()$packages)

save(packages,
     file=file.path(model.save,paste0("packages_",runName,".RData")))

### ============================================================================
### Load data
### ============================================================================

NSH.sams  <- new("FLSAMs")
NSHs      <- new("FLStocks")
NSH.ctrls <- list()

for(SMSkeyRuns in SMSkeyRunsMat){
  
  load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun
  
  NSH     <- window(NSH, end=min(SMSkeyRunsMat))
  NSH.tun$HERAS <- window(NSH.tun$HERAS, end=min(SMSkeyRunsMat))
  NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=min(SMSkeyRunsMat))
  NSH.tun$`LAI-ORSH` <- window(NSH.tun$`LAI-ORSH`, end=min(SMSkeyRunsMat))
  NSH.tun$`LAI-BUN` <- window(NSH.tun$`LAI-BUN`, end=min(SMSkeyRunsMat))
  NSH.tun$`LAI-CNS` <- window(NSH.tun$`LAI-CNS`, end=min(SMSkeyRunsMat))
  NSH.tun$`LAI-SNS` <- window(NSH.tun$`LAI-SNS`, end=min(SMSkeyRunsMat))
  NSH.tun$`IBTS-Q3` <- window(NSH.tun$`IBTS-Q3`, end=min(SMSkeyRunsMat))
  NSH.tun$`IBTS-Q1` <- window(NSH.tun$`IBTS-Q1`, end=min(SMSkeyRunsMat)+1)
  NSH.tun$IBTS0     <- window(NSH.tun$IBTS0, end=min(SMSkeyRunsMat)+1)
  
  pg <- range(NSH)[2]
  
  ### ============================================================================
  ### Construct single fleet control object
  ### ============================================================================
  
  NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)
  
  # validObject(NSH.ctrl)
  #NSH.ctrl@cor.F <- 2
  
  NSH.ctrls[[ac(SMSkeyRuns)]] <- NSH.ctrl

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
                      paste0("results_",runName,'_',SMSkeyRuns,"_sf.RData")))
  
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
  
  NSH.sams[[ac(SMSkeyRuns)]]  <- NSH.sam.opt
  NSHs[[ac(SMSkeyRuns)]]      <- NSH.opt
}

save(NSH.ctrls,
     file=file.path(model.save,
                    paste0("config_",runName,"_sf.RData")))

save(NSHs, NSH.tun, NSH.ctrls, NSH.sams,
     file=file.path(model.save,
                    paste0("results_",runName,"_sf.RData")))