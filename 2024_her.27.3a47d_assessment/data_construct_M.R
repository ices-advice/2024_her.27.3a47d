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

# General

assessment_year     <- "2024"
assessment_name     <- paste0('NSAS_HAWG', assessment_year)

### ============================================================================
### Prepare stock object for single fleet assessment
### ============================================================================

M2M1_raw     <- read.csv(file.path(data.source,"SMS_NSAS_M_raw.csv"),header=TRUE,check.names = FALSE)

### ============================================================================
### process time series
### ============================================================================

keyRuns <- unique(M2M1_raw$Source)

for(iSMS in keyRuns){
  M2M1_raw_filt <- subset(M2M1_raw,Source == iSMS & year >= 1974)
  agesUnique <- unique(M2M1_raw_filt$age)
  
  storeSmooth   <- array(NA,dim=c(length(sort(unique(M2M1_raw_filt$age))),
                                  length(sort(unique(M2M1_raw_filt$year))),
                                  3),
                         dimnames=list(age=sort(unique(M2M1_raw_filt$age)),
                                       year=sort(unique(M2M1_raw_filt$year)),
                                       Fit=c("5%","50%","95")))
  # loop on ages
  for(iAge in agesUnique){
    res         <- predict(loess((M)~year,data=subset(M2M1_raw_filt, age == iAge),span=0.5),
                           newdata=expand.grid(year=sort(unique(subset(M2M1_raw_filt,age == iAge)$year))),
                           se=T)                           
    yrs         <- sort(unique(subset(M2M1_raw_filt, age == iAge)$year))
    
    storeSmooth[ac(iAge),ac(yrs),] <- matrix(c(res$fit-1.96*res$se.fit,res$fit,res$fit+1.96*res$se.fit),nrow=length(yrs),ncol=3)
  }
  
  if(max(agesUnique) < 8){
    res <- rbind(storeSmooth[,,"50%"],storeSmooth[8,,"50%"])
    rownames(res) <- c(agesUnique,8)
  }else{
    res <- storeSmooth[,,"50%"]
  }
  
  write.csv(res,
            file=file.path(data.source,
                           paste0("M_NSAS_smoothedSpan50_notExtrapolated_sms",ac(iSMS),".csv")))
}
