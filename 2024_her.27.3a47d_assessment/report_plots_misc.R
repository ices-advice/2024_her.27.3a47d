## Prepare plots for report

## Before: results_sf.RData (model)
## After:  summary.png (report)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(ggplot2)

mkdir("report")

data.source     <- file.path(".", "bootstrap",'data')
data.save       <- file.path(".", "data")
model.save      <- file.path(".", "model")
output.save     <- file.path(".", "output")

source('utilities_plots.R')

################################################################################
### ============================================================================
### Single-fleet
### ============================================================================
################################################################################
prefix <- 'sf'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

load(file.path('../../model','NSAS_HAWG2024_sf_IBP2021.RData'))

tab.ssb <- ssb(NSH.sam)
tab.ssb$metric <- 'ssb'
tab.ssb$assessment <- 'IBP2021'
tab.rec <- rec(NSH.sam)
tab.rec$metric <- 'rec'
tab.rec$assessment <- 'IBP2021'

tab.all <- rbind(tab.ssb,tab.rec)

load(file.path('../../model','NSAS_HAWG2024_sf.RData'))

tab.ssb <- ssb(NSH.sam)
tab.ssb$metric <- 'ssb'
tab.ssb$assessment <- 'HAWG2024'
tab.rec <- rec(NSH.sam)
tab.rec$metric <- 'rec'
tab.rec$assessment <- 'HAWG2024'

tab.all <- rbind(tab.all,tab.ssb,tab.rec)

tab.all <- tab.all %>% 
            select(-c('CV','lbnd','ubnd')) %>% 
            pivot_wider(names_from = 'metric',values_from = 'value')

taf.png(paste0(prefix,"_SR_comp"))
ggplot(tab.all,aes(x=ssb,y=rec,label=year,col=assessment))+
  geom_text()+
  xlim(0,5.5e6)+
  ylim(0,7e7)
dev.off()

taf.png(paste0(prefix,"_SR_comp_recent"))
ggplot(subset(tab.all,year >= 2002),aes(x=ssb,y=rec,label=year,col=assessment))+
  geom_text()+
  xlim(1e6,3e6)+
  ylim(0,7e7)
dev.off()

### ============================================================================
### comparison with and without IBTS-Q1 age 1
### ============================================================================
NSH.sams  <- new("FLSAMs")

load(file.path('../../model','NSAS_HAWG2024_sf.RData'))
NSH.sams[['baseline']] <- NSH.sam

load(file.path('../../model','NSAS_HAWG2024_sf_no2024IBTSQ1.RData'))
NSH.sams[['no IBTS-Q1 2024']] <- NSH.sam

taf.png(paste0(prefix,"_explo_stock_trajectory_drop IBTSQ1 age 1"))
plot(NSH.sams)
dev.off()

### ============================================================================
### comparison 2024 assessment with 2023 and 2019 SMS keyruns
### ============================================================================
NSH.sams  <- new("FLSAMs")

load(file.path('../../model','NSAS_HAWG2024_sf.RData'))
NSH.sams[['baseline']] <- NSH.sam
df.m <- as.data.frame(NSH@m)
df.m$keyrun <- 'SMS2023'

load(file.path('../../model','NSAS_HAWG2024_sf_SMS2019.RData'))
NSH.sams[['SMS2019']] <- NSH.sam
temp <- as.data.frame(NSH@m)
temp$keyrun <- 'SMS2019'

df.m <- rbind(df.m,temp)

taf.png(paste0(prefix,"_explo_stock_trajectory_SMS2019-2023"))
plot(NSH.sams)
dev.off()

taf.png(paste0(prefix,"_M smooth"))
ggplot(subset(df.m,year > 1974),aes(x=year,y=data,col=keyrun))+
  geom_line()+
  facet_wrap(~age)
dev.off()

setwd('../..')

################################################################################
### ============================================================================
### multi-fleet exploration
### ============================================================================
################################################################################

load(file.path(model.save,paste0('NSAS_HAWG2024','_mf.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam
load(file.path(model.save,paste0('NSAS_HAWG2024','_mf_retro.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam

prefix <- 'mf'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

### ============================================================================
### multi-fleet comparison of age share between fleets for 1947-1996
### ============================================================================

scalA <- NSH3f@catch.n[,ac(1997:2007),,,'A']/areaSums(NSH3f@catch.n[,ac(1997:2007),,,])
scalA <- as.data.frame(scalA)
scalBD <- NSH3f@catch.n[,ac(1997:2007),,,'BD']/areaSums(NSH3f@catch.n[,ac(1997:2007),,,])
scalBD <- as.data.frame(scalBD)
scalC <- NSH3f@catch.n[,ac(1997:2007),,,'C']/areaSums(NSH3f@catch.n[,ac(1997:2007),,,])
scalC <- as.data.frame(scalC)

ageScal <- rbind(scalA,scalBD,scalC)

taf.png(paste0(prefix,"_fleet.prop"))
print(ggplot(ageScal,aes(x=year,y=data,fill=as.factor(area)))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~age))
dev.off()

### ============================================================================
### multi-fleet comparison HAWG2023
### ============================================================================

HAWG2023_ssb <- read.csv(file.path('../../bootstrap/initial/data/HAWG2023_mf_ssb.csv'))
HAWG2023_ssb$type <- 'ssb'
HAWG2023_ssb$WG <- 'HAWG2023'
HAWG2023_fbar <- read.csv(file.path('../../bootstrap/initial/data/HAWG2023_mf_fbar.csv'))
HAWG2023_fbar$type <- 'fbar'
HAWG2023_fbar$WG <- 'HAWG2023'
HAWG2023_rec <- read.csv(file.path('../../bootstrap/initial/data/HAWG2023_mf_rec.csv'))
HAWG2023_rec$type <- 'rec'
HAWG2023_rec$WG <- 'HAWG2023'

HAWG2023 <- rbind(HAWG2023_ssb,HAWG2023_fbar,HAWG2023_rec)

HAWG2024_ssb <- ssb(NSH3f.retro[['2022']])
HAWG2024_ssb$type <- 'ssb'
HAWG2024_ssb$WG <- 'HAWG2024'
HAWG2024_fbar <- fbar(NSH3f.retro[['2022']])
HAWG2024_fbar$type <- 'fbar'
HAWG2024_fbar$WG <- 'HAWG2024'
HAWG2024_rec <- rec(NSH3f.retro[['2022']])
HAWG2024_rec$type <- 'rec'
HAWG2024_rec$WG <- 'HAWG2024'

HAWG2024 <- rbind(HAWG2024_ssb,HAWG2024_fbar,HAWG2024_rec)

df.plot <- rbind(HAWG2023,HAWG2024)

taf.png(paste0(prefix,"_HAWG comp recent"))
print(ggplot(subset(df.plot,year>=2000),aes(x=year,y=value,fill=WG,col=WG))+
        geom_line()+
        geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
        facet_wrap(~type,scales='free',ncol=1))
dev.off()

taf.png(paste0(prefix,"_HAWG comp all"))
print(ggplot(df.plot,aes(x=year,y=value,fill=WG,col=WG))+
        geom_line()+
        geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
        facet_wrap(~type,scales='free',ncol=1))
dev.off()

setwd('../..')

### ============================================================================
### single-fleet comparison HAWG2023
### ============================================================================