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

data.source             <- file.path(".", "bootstrap",'data')
data.save               <- file.path(".", "data")
model.save              <- file.path(".",'model')

source('utilities_model_config.R')

### ============================================================================
### ============================================================================
### mf crazy fleet proportions
### ============================================================================
### ============================================================================

### ============================================================================
### Construct single fleet control object
### ============================================================================

runName <- 'NSAS_HAWG2024_mf_crazy'

SMSkeyRuns      <- 2023
retroFlag       <- TRUE

addM    <- 0.02

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_mf.RData")))  # NSH, NSH.tun
load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))  # NSH, NSH.tun

NSH3f@m <- NSH3f@m+addM

pg <- range(NSH3f)[2]

scalA <- NSH3f@catch.n[,ac(1997),,,'A']/areaSums(NSH3f@catch.n[,ac(1997),,,])
scalBD <- NSH3f@catch.n[,ac(1997),,,'BD']/areaSums(NSH3f@catch.n[,ac(1997),,,])
scalC <- NSH3f@catch.n[,ac(1997),,,'C']/areaSums(NSH3f@catch.n[,ac(1997),,,])

NSH3f@catch.n[,ac(1947:1996),,,"A"] <- sweep(NSH@catch.n[,ac(1947:1996)],1,scalA,"*")
NSH3f@catch.n[,ac(1947:1996),,,"BD"] <- sweep(NSH3f@catch.n[,ac(1947:1996),,,"BD"],1,scalBD,"*")
NSH3f@catch.n[,ac(1947:1996),,,"C"] <- sweep(NSH3f@catch.n[,ac(1947:1996),,,"C"],1,scalC,"*")

NSH3f@catch.wt[,ac(1947:1996),,,"A"] <- NSH3f@catch.wt[,ac(1997),,,"A"]
NSH3f@catch.wt[,ac(1947:1996),,,"BD"] <- NSH3f@catch.wt[,ac(1997),,,"BD"]
NSH3f@catch.wt[,ac(1947:1996),,,"C"] <- NSH3f@catch.wt[,ac(1997),,,"C"]
NSH3f@landings.n <- NSH3f@catch.n
NSH3f@landings.wt <- NSH3f@catch.wt

NSH3f@name     <- runName

NSH3f@catch        <- computeCatch(NSH3f)

NSH3.ctrl <- config_mf_IBPNSherring2021_alt(NSH3f,NSH.tun,pg)
NSH3.ctrl@residuals <- TRUE

#NSHs3$residual@catch.n[ac(0),ac(2022),,,'A'] <- 0

### ============================================================================
### Run multi fleet model
### ============================================================================

# Run model
NSH3f.sam   <- FLSAM(NSH3f,
                     NSH.tun,
                     NSH3.ctrl)

# sam.fit <- FLSAM(NSHs3,
#                  NSH.tun,
#                  NSH3.ctrl,
#                  starting.values = NSH3f.sam.stk0,return.fit=TRUE)

save(NSH3f,
     NSH.tun,
     NSH3.ctrl,
     NSH3f.sam,
     file=file.path(model.save,
                    paste0(runName,".RData")))

### ============================================================================
### run multi fleet retro
### ============================================================================

if(retroFlag){
  # retro
  n.retro.years <- 10  # Number of years for which to run the retrospective
  NSH3.ctrl@residuals <- FALSE
  NSH3f.retro <- retro(NSH3f, NSH.tun, NSH3.ctrl, retro=n.retro.years)
  
  save(NSH3f.retro,
       file=file.path(model.save,paste0(runName,'_retro.RData')))
}

### ============================================================================
### Run multi fleet retro
### ============================================================================
load(file.path(model.save,paste0('NSAS_HAWG2024','_sf.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam

NSH.sams  <- new("FLSAMs")
NSH.sams[['sf']] <- NSH.sam
NSH.sams[['mf']] <- NSH3f.sam

out_assessment <- rbind(cbind(ssb(NSH.sams),type='ssb'),
                        cbind(rec(NSH.sams),type='rec'),
                        cbind(fbar(NSH.sams),type='fbar'))

windows()
ggplot(subset(out_assessment,year>=2000),aes(x=year,y=value,fill=name,col=name))+
  geom_line()+
  geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
  facet_wrap(~type,scales='free',ncol=1)

windows()
ggplot(out_assessment,aes(x=year,y=value,fill=name,col=name))+
  geom_line()+
  geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
  facet_wrap(~type,scales='free',ncol=1)
