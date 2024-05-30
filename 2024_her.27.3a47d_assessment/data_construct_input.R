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
SMSkeyRunsMat     <- c('2010','2013','2016','2019','2023')
assessment_year   <- "2024"
assessment_name   <- paste0('NSAS_HAWG', assessment_year)

for(SMSkeyRuns in SMSkeyRunsMat){
  
  print(paste("SMS key run:", SMSkeyRuns ))
  
  ### ============================================================================
  ### Prepare stock object for single fleet assessment
  ### ============================================================================
  
  ## Load object
  #NSH <- readFLStock("index.txt", no.discards=TRUE)
  NSH <- readFLStock(file.path(data.source,"index.txt"),no.discards=TRUE,quiet=FALSE)
  
  #Catch is calculated from: catch.wt * catch.n, however, the reported landings are
  #normally different (due to SoP corrections). Hence we overwrite the calculate landings
  # are we not using catches data then?
  NSH@catch           <- NSH@landings
  units(NSH)[1:17]    <- as.list(c(rep(c("tonnes","thousands","kg"),4), 
                                   rep("NA",2),"f",rep("NA",2)))
  
  #Set object details
  NSH@name                              <- paste0(assessment_name,'_sf')
  range(NSH)[c("minfbar","maxfbar")]    <- c(2,6)
  NSH                                   <- setPlusGroup(NSH,NSH@range["max"])
  
  #Historical data is only provided for ages 0-8 prior to 1960 (rather than 0-9 after)
  #We therefore need to fill in the last ages by applying the following assumptions
  #  - weight in the catch and weight in the stock at age 9 is the same as 
  #    in period 1960-1983 (constant in both cases)
  #  - catch at age reported as 8+ is split evenly between age 8 and 9
  #  - natural mortality at age 9 is 0.1 (same as age 8)
  #  - proportion mature at age 9 is 1.0 (same as age 8)
  #  - harvest.spwn and m.spwn are the same as elsewhere
  hist.yrs                      <- as.character(1947:1959)
  NSH@catch.wt                  <- NSH@landings.wt #Automatic population of catch.wt introduces NAS
  NSH@catch.wt["9",hist.yrs]    <- 0.271
  NSH@landings.wt["9",hist.yrs] <- 0.271
  NSH@catch.n["9",hist.yrs]     <- NSH@catch.n["8",hist.yrs]/2
  NSH@catch.n["8",hist.yrs]     <- NSH@catch.n["9",hist.yrs]
  NSH@landings.n["9",hist.yrs]  <- NSH@landings.n["8",hist.yrs]/2
  NSH@landings.n["8",hist.yrs]  <- NSH@landings.n["9",hist.yrs]
  NSH@stock.wt["9",hist.yrs]    <- 0.312
  NSH@m["9",hist.yrs]           <- 0.1
  NSH@mat["9",hist.yrs]         <- 1
  
  #No catches of age 9 in 1977 so stock.wt does not get filled there.
  #Hence, we copy the stock weight for that age from the previous year.
  #Note that because we use a fixed stock.wt prior to 1983, there is no
  #need to use averaging or anything fancier
  NSH@stock.wt["9","1977"]              <- NSH@stock.wt["9","1976"]
  
  #Use a running mean(y-2,y-1,y) of input wests (i.e. west_raw) to calculate west
  NSH@stock.wt[,3:dim(NSH@stock.wt)[2]] <- (NSH@stock.wt[,3:(dim(NSH@stock.wt)[2]-0)] +
                                              NSH@stock.wt[,2:(dim(NSH@stock.wt)[2]-1)] +
                                              NSH@stock.wt[,1:(dim(NSH@stock.wt)[2]-2)]) / 3
  
  ### ============================================================================
  ### Prepare natural mortality estimates
  ### ============================================================================
  
  ## Read in estimates from external file
  M2           <- read.csv(file.path(data.source,
                                     paste0("M_NSAS_smoothedSpan50_notExtrapolated_sms",ac(SMSkeyRuns),".csv")))
  colnames(M2)  <- sub("X","",colnames(M2))
  rownames(M2)  <- M2[,1]
  M2            <- M2[,-1]# Trim off first column as it contains 'ages'
  M2            <- M2[,apply(M2,2,function(x){all(is.na(x))==F})] # keep only years with data
  
  
  #Extract key data from default assessment
  NSHM2           <- NSH
  NSHM2@m[]       <- NA
  yrs             <- dimnames(NSHM2@m)$year
  yrs             <- yrs[which(yrs %in% colnames(M2))]
  NSHM2@m[rownames(M2),yrs][] <- as.matrix(M2)
  
  #- M extension to previous and future years (relative to vector given by WGSAM). 
  # One uses a 5 year running average for that
  
  # year discrepency between assessment and M from WGSAM
  extryrs         <- dimnames(NSHM2@m)$year[which(!dimnames(NSHM2@m)$year %in% yrs)]
  # year after those avaible for M from WGSAM
  extryrsfw       <- extryrs[which(extryrs > max(an(yrs)))]
  # years prior to those available for M from WGSAM
  extryrsbw       <- extryrs[which(extryrs <= max(an(yrs)))]
  ages            <- dimnames(NSHM2@m)$age
  extrags         <- names(which(apply(M2,1,function(x){all(is.na(x))==T})==T))
  yrAver          <- 5
  for(iYr in as.numeric(rev(extryrs))){
    for(iAge in ages[!ages%in%extrags]){
      if(iYr %in% extryrsbw){
        NSHM2@m[ac(iAge),ac(iYr)] <- yearMeans(NSHM2@m[ac(iAge),ac((iYr+1):(iYr+yrAver)),],na.rm=T)
      }
      if(iYr %in% extryrsfw){
        NSHM2@m[ac(iAge),ac(iYr)] <- yearMeans(NSHM2@m[ac(iAge),ac((iYr-1):(iYr-yrAver)),],na.rm=T)
      }
    }
  }
  if(length(extrags)>0){
    for(iAge in extrags)
      NSHM2@m[ac(iAge),]          <- NSHM2@m[ac(as.numeric(min(sort(extrags)))-1),]
  }
  #Write new M values into the original stock object
  #addM      <- 0.11 #M profiling based on 2018 benchmark meeting
  #addM      <- 0
  NSH@m     <- NSHM2@m
  
  ### ============================================================================
  ### Prepare index object for assessment
  ### ============================================================================
  
  ## Load and modify all numbers at age data
  NSH.tun   <- readFLIndices(file.path(data.source,"fleet.txt"))
  surveyLAI <- read.table(file.path(data.source,"lai.txt"), stringsAsFactors=FALSE, header=TRUE)
  
  NSH.tun   <- lapply(NSH.tun,function(x) {x@type <- "number"; return(x)})
  NSH.tun[["IBTS0"]]@range["plusgroup"] <- NA
  
  # subset LAI index
  ORSH          <- subset(surveyLAI,Area == "Or/Sh")
  CNS           <- subset(surveyLAI,Area == "CNS")
  BUN           <- subset(surveyLAI,Area == "Buchan")
  SNS           <- subset(surveyLAI,Area == "SNS")
  
  ORSH                <- formatLAI(ORSH,1972,range(NSH)["maxyear"])
  CNS                 <- formatLAI(CNS,1972,range(NSH)["maxyear"])
  BUN                 <- formatLAI(BUN,1972,range(NSH)["maxyear"])
  SNS                 <- formatLAI(SNS,1972,range(NSH)["maxyear"])
  
  FLORSH              <- FLIndex(index=FLQuant(t(ORSH),dimnames=list(age=colnames(ORSH),year=rownames(ORSH),unit="ORSH",season="all",area="unique",iter="1")))
  FLCNS               <- FLIndex(index=FLQuant(t(CNS),dimnames=list(age=colnames(CNS),year=rownames(CNS),unit="CNS",season="all",area="unique",iter="1")))
  FLBUN               <- FLIndex(index=FLQuant(t(BUN),dimnames=list(age=colnames(BUN),year=rownames(BUN),unit="BUN",season="all",area="unique",iter="1")))
  FLSNS               <- FLIndex(index=FLQuant(t(SNS),dimnames=list(age=colnames(SNS),year=rownames(SNS),unit="SNS",season="all",area="unique",iter="1")))
  range(FLORSH)[6:7]  <- range(FLCNS)[6:7] <- range(FLBUN)[6:7] <- range(FLSNS)[6:7] <- c(0.67,0.67)
  name(FLORSH)        <- "LAI-ORSH"
  name(FLCNS)         <- "LAI-CNS"
  name(FLBUN)         <- "LAI-BUN"
  name(FLSNS)         <- "LAI-SNS"
  type(FLORSH)        <- type(FLCNS) <- type(FLBUN) <- type(FLSNS) <- "partial"
  FLORSH@index@.Data[which(is.na(FLORSH@index))] <- -1
  FLCNS@index@.Data[which(is.na(FLCNS@index))] <- -1
  FLBUN@index@.Data[which(is.na(FLBUN@index))] <- -1
  FLSNS@index@.Data[which(is.na(FLSNS@index))] <- -1
  NSH.tun[(length(NSH.tun)+1):(length(NSH.tun)+4)] <- c(FLORSH,FLBUN,FLCNS,FLSNS)
  names(NSH.tun)[(length(NSH.tun)-3):(length(NSH.tun))] <- paste("LAI",c("ORSH","BUN","CNS","SNS"),sep="-")
  
  ### ============================================================================
  ### Apply plus group to all data sets
  ### ============================================================================
  
  pg <- max(an(rownames(M2)))
  if(pg > 8) pg <- 8
  #pg <- 8
  
  #- This function already changes the stock and landings.wts correctly
  NSH <- setPlusGroup(NSH,pg)
  
  NSH.tun[["HERAS"]]@index[ac(pg),]     <- quantSums(NSH.tun[["HERAS"]]@index[ac(pg:dims(NSH.tun[["HERAS"]]@index)$max),])
  NSH.tun[["HERAS"]]                    <- trim(NSH.tun[["HERAS"]],age=dims(NSH.tun[["HERAS"]]@index)$min:pg)
  NSH.tun[["HERAS"]]@range["plusgroup"] <- pg
  
  NSH.tun[["IBTS-Q3"]] <- trim(NSH.tun[["IBTS-Q3"]],age=0:5)
  NSH.tun[["IBTS-Q1"]] <- trim(NSH.tun[["IBTS-Q1"]],age=1)
  NSH.tun[["HERAS"]] <- trim(NSH.tun[["HERAS"]],age=1:pg)
  
  
  ### ============================================================================
  ### Closure data deletion
  ### ============================================================================
  
  ## We don't believe the closure catch data, so put it to NA
  NSH@catch.n[,ac(1978:1979)] <- NA
  
  ### ============================================================================
  ### Now the multifleet data
  ### ============================================================================

  # read the catch data
  caaF  <- read.table(file.path(data.source, "canum_mf.txt"),stringsAsFactors=F,header=T)
  caA   <- matrix(subset(caaF,Fleet=="A")[,"numbers"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))*1000
  caB   <- matrix(subset(caaF,Fleet=="B")[,"numbers"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))*1000
  caC   <- matrix(subset(caaF,Fleet=="C")[,"numbers"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))*1000
  caD   <- matrix(subset(caaF,Fleet=="D")[,"numbers"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))*1000
  cwA   <- matrix(subset(caaF,Fleet=="A")[,"weight"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))
  cwB   <- matrix(subset(caaF,Fleet=="B")[,"weight"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))
  cwC   <- matrix(subset(caaF,Fleet=="C")[,"weight"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))
  cwD   <- matrix(subset(caaF,Fleet=="D")[,"weight"],nrow=length(0:9),ncol=length(1997:range(NSH)["maxyear"]),dimnames=list(age=0:9,year=1997:range(NSH)["maxyear"]))

  #- Functions to set the plusgroups
  setMatrixPlusGroupN <- function(x,pg){
    idxpg <- which(an(rownames(x)) == pg)
    x[idxpg,] <- colSums(x[idxpg:nrow(x),],na.rm=T)
    x <- x[1:idxpg,]
    return(x)}
  
  setMatrixPlusGroupWt <- function(x,y,pg){
    idxpg <- which(an(rownames(x)) == pg)
    y[idxpg,] <- colSums(x[idxpg:(nrow(x)),]*y[idxpg:nrow(y),],na.rm=T) / colSums(x[idxpg:nrow(x),],na.rm=T)
    y[is.nan(y)] <- 0
    y         <- y[1:idxpg,]
    return(y)}
  
  NSH3f <- FLCore::expand(NSH,area=c("A","BD","C"))
  NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"A"]   <- setMatrixPlusGroupN(caA,pg)
  NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"BD"]   <- setMatrixPlusGroupN(caB+caD,pg)
  NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"C"]   <- setMatrixPlusGroupN(caC,pg)
  NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"A"]  <- setMatrixPlusGroupWt(caA,cwA,pg)
  NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"BD"]  <- setMatrixPlusGroupWt(caB+caD,(cwB*caB + cwD*caD)/(caB+caD),pg)
  NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"C"]  <- setMatrixPlusGroupWt(caC,cwC,pg)
  
  scalA <- NSH3f@catch.n[,ac(1997),,,'A']/areaSums(NSH3f@catch.n[,ac(1997),,,])
  scalBD <- NSH3f@catch.n[,ac(1997),,,'BD']/areaSums(NSH3f@catch.n[,ac(1997),,,])
  scalC <- NSH3f@catch.n[,ac(1997),,,'C']/areaSums(NSH3f@catch.n[,ac(1997),,,])
  
  NSH3f@catch.n[,ac(1947:1996),,,"A"] <- sweep(NSH@catch.n[,ac(1947:1996)],1,scalA,"*")
  NSH3f@catch.n[,ac(1947:1996),,,"BD"] <- sweep(NSH@catch.n[,ac(1947:1996)],1,scalBD,"*")
  NSH3f@catch.n[,ac(1947:1996),,,"C"] <- sweep(NSH@catch.n[,ac(1947:1996)],1,scalC,"*")
  
  NSH3f@catch.wt[,ac(1947:1996),,,"A"] <- NSH3f@catch.wt[,ac(1997),,,"A"]
  NSH3f@catch.wt[,ac(1947:1996),,,"BD"] <- NSH3f@catch.wt[,ac(1997),,,"BD"]
  NSH3f@catch.wt[,ac(1947:1996),,,"C"] <- NSH3f@catch.wt[,ac(1997),,,"C"]
  NSH3f@landings.n <- NSH3f@catch.n
  NSH3f@landings.wt <- NSH3f@catch.wt
  
  NSH3f@name     <- paste0(assessment_name,'_mf')
  
  NSH3f@catch        <- computeCatch(NSH3f)
  
  NSH3f@catch.n[,ac(1978:1979)] <- NA

  range(NSH3f)[c("minfbar","maxfbar")]     <- c(2,6)
  
  ### ============================================================================
  ### Save objects
  ### ============================================================================

  save(NSH, NSH.tun, 
       file=file.path(data.save,paste0("data_sms",SMSkeyRuns,"_sf.RData")))
  
  save(NSH3f, NSH.tun, 
       file=file.path(data.save,paste0("data_sms",SMSkeyRuns,"_mf.RData")))
}

