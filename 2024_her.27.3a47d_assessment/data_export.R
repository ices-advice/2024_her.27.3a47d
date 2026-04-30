## Preprocess data, write TAF data tables
## Part 2: Export TAF tables

## Before: canum.txt, fleet.txt (bootstrap/data), data_sf.RData (data)
## After:  catage.csv, catage_full.csv, maturity.csv, natmort.csv, propf.csv,
##         propm.csv, survey_heras.csv, survey_heras_full.csv,
##         survey_ibts_0.csv, survey_ibts_q1.csv, survey_ibts_q3.csv,
##         survey_lai_bun.csv, survey_lai_cns.csv, survey_lai_orsh.csv,
##         survey_lai_sns.csv, wcatch.csv, wstock.csv (data)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)

SMSkeyRuns  <- 2019

data.save       <- file.path(".", "data")
data.source     <- file.path(".", "bootstrap",'data')

### ============================================================================
### Read data
### ============================================================================

## Full tables
catage.full <- flr2taf(readVPAFile(file.path(data.source,'canum.txt')))
survey.full <- readFLIndices(file.path(data.source,'fleet.txt'))
survey.heras.full <- flr2taf(survey.full$HERAS@index)

## Access tables from FLR objects
load(file.path(file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData"))))

### ============================================================================
### Export tables
### ============================================================================

setwd("data")
write.taf(catage.full)
write.taf(survey.heras.full)

## Stock tables
write.taf(plus(flr2taf(NSH@catch.n)), "catage.csv")
write.taf(flr2taf(NSH@mat), "maturity.csv")
write.taf(flr2taf(NSH@m), "natmort.csv")
write.taf(flr2taf(NSH@harvest.spwn), "propf.csv")
write.taf(flr2taf(NSH@m.spwn), "propm.csv")
write.taf(flr2taf(NSH@catch.wt), "wcatch.csv")
write.taf(flr2taf(NSH@stock.wt), "wstock.csv")

## Survey tables
write.taf(plus(flr2taf(NSH.tun[[1]]@index)), "survey_heras.csv")
write.taf(flr2taf(NSH.tun[[2]]@index), "survey_ibts_q1.csv")
write.taf(flr2taf(NSH.tun[[3]]@index), "survey_ibts_0.csv")
write.taf(flr2taf(NSH.tun[[4]]@index), "survey_ibts_q3.csv")
write.taf(flr2taf(NSH.tun[[5]]@index), "survey_lai_orsh.csv")
write.taf(flr2taf(NSH.tun[[6]]@index), "survey_lai_bun.csv")
write.taf(flr2taf(NSH.tun[[7]]@index), "survey_lai_cns.csv")
write.taf(flr2taf(NSH.tun[[8]]@index), "survey_lai_sns.csv")

setwd("..")
