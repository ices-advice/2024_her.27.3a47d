## Run analysis, write model results

## Before: input.RData (data)
## After:  forecast.RData, transcript.txt (model)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(methods)
library(minpack.lm)  # nls.lm, in fleet.harvest
library(tidyverse)

mkdir("model")

sourceDir("bootstrap/initial/software/utilities")

data.source <- file.path('.','bootstrap','data')
report.dir <- file.path('.','report')

assessment_year     <- '2024'

## -----------------------------------------------------------------------------
## 1  Read forecast inputt
## -----------------------------------------------------------------------------

load(file.path('.','data',paste0('input_',assessment_year,'_nego.RData')))

DtY   <- ac(range(stf)["maxyear"]-3) #Data year
ImY   <- ac(an(DtY)+1)             #Intermediate year
FcY   <- ac(an(DtY)+2)             #Forecast year
CtY   <- ac(an(DtY)+3)             #Continuation year
CtY1  <- ac(an(DtY)+4)
FuY   <- c(ImY,FcY,CtY)            #Future years

dsc         <- "North Sea Herring"
nam         <- "NSH"
dms         <- dimnames(stf)

#--------------------------------------------
# Advice scenario
#--------------------------------------------

if("Headline advice corrected" %in% stf.options){
  
  caseName <- "Headline advice corrected"
  
  TACS.temp <- TACS
  TACS.temp[,ac(2025:2026),'D'] <- 0.1
  TACS.temp[,ac(2025:2026),'C'] <- 0.1
  
  res <- fmsyAR_fun_notransfer(  stf,
                                 FuY,
                                 TACS.temp,
                                 RECS,
                                 referencePoints,
                                 TAC_var,
                                 f01,
                                 f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                  <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]  <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                    <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                    <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]     <- res$stf@catch[,FcY]
  stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]     <- res$ssb.FcY
  stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]     <- res$ssb.CtY
}

#--------------------------------------------
# Rollover (2024 B fleet level)
#--------------------------------------------

if('C-TAC rule_B rollover_0% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B rollover_0% transfer'
  
  res <- fmsyAR_nego_TACrule_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        7716,
                                        c(1,1),
                                        c(0,0),
                                        f01,
                                        f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

if('C-TAC rule_B rollover_100% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B rollover_100% transfer'
  
  res <- fmsyAR_nego_TACrule_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        7716,
                                        c(1,1),
                                        c(1,1),
                                        f01,
                                        f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}


if('C-TAC rule_B rollover_NSAS share transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B rollover_NSAS share transfer'
  
  res <- fmsyAR_nego_TACruleTrans_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        7716,
                                        c(1,1),
                                        c(TAC_var$CtransferFcY,TAC_var$DtransferFcY),
                                        f01,
                                        f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

#--------------------------------------------
# B fleet decrease with A flee
#--------------------------------------------

if('C-TAC rule_B decrease_0% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B decrease_0% transfer'
  
  tar <- 7716-7716*0.228
  
  res <- fmsyAR_nego_TACrule_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        tar,
                                        c(1,1),
                                        c(0,0),
                                        f01,
                                        f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

if('C-TAC rule_B decrease_100% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B decrease_100% transfer'
  
  tar <- 7716-7716*0.228
  
  res <- fmsyAR_nego_TACrule_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        tar,
                                        c(1,1),
                                        c(1,1),
                                        f01,
                                        f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

#--------------------------------------------
# B fleet level as proportion of A fleet
#--------------------------------------------

if('C-TAC rule_B prop_0% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B prop_0% transfer'
  
  catchBtarProp <- 7716/510323
  
  res <- fmsyAR_nego_TACrule_catchBtarProp( stf,
                                            FuY,
                                            TACS,
                                            RECS,
                                            referencePoints,
                                            TAC_var,
                                            catchBtarProp,
                                            c(0,0),
                                            f01,
                                            f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}


if('C-TAC rule_B prop_100% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_B prop_100% transfer'
  
  catchBtarProp <- 7716/510323
  
  res <- fmsyAR_nego_TACrule_catchBtarProp( stf,
                                            FuY,
                                            TACS,
                                            RECS,
                                            referencePoints,
                                            TAC_var,
                                            catchBtarProp,
                                            c(1,1),
                                            f01,
                                            f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

#--------------------------------------------
# scale B fleet to achieve F01=0.05
#--------------------------------------------

if('C-TAC rule_F01=0.05_0% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_F01=0.05_0% transfer'
  
  res <- fmsyAR_nego_TACrule_F01tar( stf,
                                     FuY,
                                     TACS,
                                     RECS,
                                     referencePoints,
                                     TAC_var,
                                     0.05,
                                     c(0,0),
                                     f01,
                                     f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}


if('C-TAC rule_F01=0.05_100% transfer' %in% stf.options){
  
  caseName <- 'C-TAC rule_F01=0.05_100% transfer'
  
  res <- fmsyAR_nego_TACrule_F01tar( stf,
                                     FuY,
                                     TACS,
                                     RECS,
                                     referencePoints,
                                     TAC_var,
                                     0.05,
                                     c(1,1),
                                     f01,
                                     f26)
  
  # update stf table
  stf.table[caseName,"Fbar 2-6 A",]                                   <- quantMeans(res$stf@harvest[f26,FcY,"A"])
  stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),]   <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
  stf.table[caseName,"Fbar 1-3 C",]                                   <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
  stf.table[caseName,"Fbar 2-6",]                                     <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
  stf.table[caseName,"Fbar 0-1",]                                     <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
  stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]      <- res$stf@catch[,FcY]
  
  stf.table[caseName,'TAC A+Ctrans',]      <- res$stf@catch[,FcY,'A']
  stf.table[caseName,'TAC A only',]      <- res$TAC_opti[,FcY,'A']-res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-NS',]      <- res$TAC_opti[,FcY,'C_NS']
  stf.table[caseName,'TAC C-3a',]      <- res$TAC_opti[,FcY,'C_3a']
  
  stf.table[caseName,'TAC B+Dtrans',]      <- res$stf@catch[,FcY,'B']
  stf.table[caseName,'TAC B only',]      <- res$stf@catch[,FcY,'B']-res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-NS',]      <- res$TAC_opti[,FcY,'D_NS']
  stf.table[caseName,'TAC D-3a',]      <- res$TAC_opti[,FcY,'D_3a']
  
  stf.table[caseName,'prop_C_NS',]      <- res$prop_C_NS
  
  stf.table[caseName,"SSB_FcY",]      <- res$ssb.FcY
  stf.table[caseName,"SSB_CtY",]      <- res$ssb.CtY
}

## -----------------------------------------------------------------------------
## 3  saving outputs
## -----------------------------------------------------------------------------

# saving objects
save(stf, stf.table, 
     file=file.path('.','model',paste0('forecast_',assessment_year,'_nego.RData')))

# output tables
forecast <- xtab2taf(stf.table[,,"50%"])
# names(forecast) <- c("Basis", "Fbar26A", "Fbar01B", "Fbar13C", "Fbar01D",
#                      "Fbar26", "Fbar01", "CatchA", "TAC A","CatchB","TAC B", "CatchC","TAC C", "CatchD","TAC D",
#                      "SSB1", "SSB2")

write.taf(forecast, dir="output",file = paste0('forecast_',assessment_year,'_nego.csv'))

## -----------------------------------------------------------------------------
## scanning uptake levels
## -----------------------------------------------------------------------------

uptakeVec <- c(0.5,0.6,0.7,0.8,0.9,1)
df.out <- data.frame(matrix(ncol = 3, nrow = length(uptakeVec)))
colnames(df.out) <- c('uptake','BD_catchNS','F01')

for(uptakeIdx in 1:length(uptakeVec)){
  print(uptakeVec[uptakeIdx])
  res <- fmsyAR_nego_TACrule_catchBtar( stf,
                                        FuY,
                                        TACS,
                                        RECS,
                                        referencePoints,
                                        TAC_var,
                                        7716,
                                        c(uptakeVec[uptakeIdx],uptakeVec[uptakeIdx]),
                                        c(1,1),
                                        f01,
                                        f26)
  
  df.out$uptake[uptakeIdx] <- uptakeVec[uptakeIdx]
  df.out$BD_catchNS[uptakeIdx] <- res$stf@catch[,FcY,'B']+res$stf@catch[,FcY,'D']
  df.out$F01[uptakeIdx] <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
}

p <- ggplot(df.out, aes(x=BD_catchNS,y=F01))+
      theme_bw()+
      ylab('F01')+
      geom_hline(yintercept = 0.05,col='red')+
      xlab('Catch BD (t)')+
      geom_line()

ggsave(file.path(report.dir,'scan BD F01 vs catches.png'),
       p,
       width = 170,
       height = 100,
       units = c("mm"),
       dpi = 300)

coeff <- mean(df.out$BD_catchNS/df.out$F01)

p <- ggplot(df.out, aes(x=uptake)) +
      theme_bw()+
      geom_line( aes(y=F01),color='blue') + 
      geom_line( aes(y=BD_catchNS / coeff),color='green')+ # Divide by 10 to get the same range than the temperature
      geom_hline(yintercept = 0.05,col='red')+
      scale_y_continuous(name = "F01",
                         sec.axis = sec_axis(~.*coeff, name="Catch (t)"))+
      theme(axis.title.y = element_text(color = 'blue'),
            axis.title.y.right = element_text(color = 'green'))

ggsave(file.path(report.dir,'scan uptake BD.png'),
       p,
       width = 170,
       height = 100,
       units = c("mm"),
       dpi = 300)
