## Preprocess data, write TAF data tables
## Part 1: Construct FLR objects

## Before: canum.txt, caton.txt, fleet.txt, fprop.txt, index.txt, lai.txt,
##         matprop.txt, mprop.txt,
##         Smoothed_span50_M_NotExtrapolated_NSASSMS2016.csv, weca.txt,
##         west_raw.txt (bootstrap/data)
##         canum
## After:  data_sf.RData (data)

### ============================================================================
### init
### ============================================================================

rm(list=ls())

library(icesTAF)
library(tidyr)

taf.library(FLCore)

# sessionInfo()
# taf.session()

library(methods)

source("utilities_data.R")

mkdir("data")

data.source   <- file.path(".","bootstrap", "data")
data.save     <- file.path(".", "data")

# General
SMSkeyRuns     <- 2023
assessment_year   <- "2024"
assessment_name   <- paste0('NSAS_HAWG', assessment_year)

load(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_mf.RData")))  # NSH, NSH.tun

### ============================================================================
### scenario
### ============================================================================

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

NSH3f@name     <- paste0(assessment_name,'_mf_alt')

NSH3f@catch        <- computeCatch(NSH3f)

NSH3f@catch.n[,ac(1978:1979)] <- NA

save(NSH3f, NSH.tun, 
     file=file.path(data.save,paste0("data_sms",SMSkeyRuns,"_mf_alt.RData")))
