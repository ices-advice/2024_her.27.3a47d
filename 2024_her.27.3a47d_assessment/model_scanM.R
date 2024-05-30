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

runName     <- 'NSAS_HAWG2024_sf_scanM'

SMSkeyRun <- 2023

source('utilities_model.R')
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

load(file.path(data.save,paste0("data_sms",SMSkeyRun,"_sf.RData")))

pg <- range(NSH)[2]

NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,pg)

save(NSH.ctrl, 
     file=file.path(model.save,
                    paste0("config_",runName,".RData")))

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

# store optimal result
df.nlogl       <- as.data.frame(unlist(lapply(res$NSH.sams,nlogl)))
colnames(df.nlogl) <- 'nlogl'
df.nlogl$addM  <- rownames(df.nlogl)

idxOpt <- which(min(df.nlogl$nlogl) == df.nlogl$nlogl)

df.nlogl$addM[idxOpt]

NSH.sam.opt   <- res$NSH.sams[[ac(df.nlogl$addM[idxOpt])]]
M             <- res$NSH.m[[ac(df.nlogl$addM[idxOpt])]]

NSH.opt           <- NSH
NSH.opt@stock.n   <- NSH.sam.opt@stock.n[,ac(range(NSH.opt)["minyear"]:range(NSH.opt)["maxyear"])]
NSH.opt@harvest   <- NSH.sam.opt@harvest[,ac(range(NSH.opt)["minyear"]:range(NSH.opt)["maxyear"])]
NSH.opt@m         <- M

save(NSH.opt,NSH.sam.opt,NSH, NSH.tun, NSH.ctrl,res,
     file=file.path(model.save,
                    paste0(runName,".RData")))