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

assessment_year     <- '2024'

TAC_CD_advice <- T

for(TAC_CD_advice in c(TRUE,FALSE)){

  stf.res <- tibble()
  
  ## -----------------------------------------------------------------------------
  ## 1  Read forecast inputt
  ## -----------------------------------------------------------------------------
  
  if(TAC_CD_advice){
    load(file.path('.','data',paste0('input_',assessment_year,'.RData')))
  }else{
    load(file.path('.','data',paste0('input_',assessment_year,'_CDsq.RData')))
  }
  
  DtY   <- ac(range(stf)["maxyear"]-3) #Data year
  ImY   <- ac(an(DtY)+1)             #Intermediate year
  FcY   <- ac(an(DtY)+2)             #Forecast year
  CtY   <- ac(an(DtY)+3)             #Continuation year
  CtY1  <- ac(an(DtY)+4)
  FuY   <- c(ImY,FcY,CtY)            #Future years
  
  dsc         <- "North Sea Herring"
  nam         <- "NSH"
  dms         <- dimnames(stf)
  
  #-------------------------------------------------------------------------------
  # 7a) Fmsy Advice rule transfer
  #-------------------------------------------------------------------------------
  
  if("fmsyAR_TACrule_transfer" %in% stf.options){
    caseName <- "fmsyAR_TACrule_transfer"
    
    res <- fmsyAR_fun_TACrule_transfer( stf,
                                        FuY,
                                        TACS,
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
    #stf.table[caseName,"C transfer",]                                  <- res$prop_C_NS
    
    stf[,c(FcY,CtY)] <- res$stf[,c(FcY,CtY)]
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  if("fmsyAR_TACrule_notransfer" %in% stf.options){
    caseName <- "fmsyAR_TACrule_notransfer"
    
    res <- fmsyAR_fun_TACrule_notransfer( stf,
                                      FuY,
                                      TACS,
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
    #stf.table[caseName,"C transfer",]                                  <- res$prop_C_NS
    
    stf[,c(FcY,CtY)] <- res$stf[,c(FcY,CtY)]
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  if("fmsyAR" %in% stf.options){
    
    caseName <- "fmsyAR"
    
    res <- fmsyAR_fun( stf,
                       FuY,
                       TACS,
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
    #stf.table[caseName,"C transfer",]                                  <- res$prop_C_NS
    
    stf[,c(FcY,CtY)] <- res$stf[,c(FcY,CtY)]
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7a) Fmsy Advice rule with C fleet transfer and 0.05 target for F01
  #-------------------------------------------------------------------------------
  
  if("fmsyAR_Btarget" %in% stf.options){
    
    caseName <- "fmsyAR_Btarget"
    
    res <- fmsyAR_fun_Btarget( stf,
                               FuY,
                               TACS,
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
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7a) Fmsy Advice rule no transfer
  #-------------------------------------------------------------------------------
  #fmsyAR_fun_no_transfer.r
  if("fmsyAR_notransfer" %in% stf.options){
    
    caseName <- "fmsyAR_notransfer"
    
    res <- fmsyAR_fun_notransfer(  stf,
                                    FuY,
                                    TACS,
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
    #stf.table[caseName,"C transfer",]                                  <- res$prop_C_NS
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7a) Fmsy Advice rule no transfer with F01=0.05 as target
  #-------------------------------------------------------------------------------
  #fmsyAR_fun_no_transfer.r
  if("fmsyAR_Btarget_notransfer" %in% stf.options){
    
    caseName <- "fmsyAR_Btarget_notransfer"
    
    res <- fmsyAR_fun_Btarget_notransfer(  stf,
                                            FuY,
                                            TACS,
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
    
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  
  #-------------------------------------------------------------------------------
  # 7c) No fishing
  #-------------------------------------------------------------------------------
  
  if("nf" %in% stf.options){
    
    caseName <- "nf"
    
    res <- nf_fun(stf,
                  RECS,
                  FuY)
    
    # update stf table
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7d) 15% increase in TAC for the A-fleet
  #-------------------------------------------------------------------------------
  
  if("+15%" %in% stf.options){
    
    caseName <- "+15%"
    
    res <- TAC_scaling_fun( stf,
                            FuY,
                            TACS,
                            CATCH,
                            RECS,
                            TAC_var,
                            1.15)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7e) 15% reduction in TAC for the A-fleet
  #-------------------------------------------------------------------------------
  
  if("-15%" %in% stf.options){
    
    caseName <- "-15%"
    
    res <- TAC_scaling_fun( stf,
                            FuY,
                            TACS,
                            CATCH,
                            RECS,
                            TAC_var,
                            0.85)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7f) Same TAC for A-fleet as last year
  #
  # for 2015: use estimated B-fleet TAC from mp for the FcY and the NSAS 
  # proportion of the advised C-fleet catch in the FcY
  #-------------------------------------------------------------------------------
  
  if("tacro" %in% stf.options){
    
    caseName <- "tacro"
    
    res <- TAC_scaling_fun( stf,
                            FuY,
                            TACS,
                            CATCH,
                            RECS,
                            TAC_var,
                            1)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7g) Fmsy option
  # Note: it is important to have a catch assumption for the B fleet in FcY and CtY
  #-------------------------------------------------------------------------------
  
  if("fmsy" %in% stf.options){
    
    caseName <- "fmsy"
    
    res <- F_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Fmsy,
                          f01,
                          f26)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7h) Fpa option
  #-------------------------------------------------------------------------------
  
  if("fpa" %in% stf.options){
    
    caseName <- "fpa"
    
    res <- F_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Fpa,
                          f01,
                          f26)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7i) Flim option
  #-------------------------------------------------------------------------------
  
  if("flim" %in% stf.options){
    
    caseName <- "flim"
    
    res <- F_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Flim,
                          f01,
                          f26)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7j) Fsq option
  #-------------------------------------------------------------------------------
  
  if("fsq" %in% stf.options){
    
    caseName <- "fsq"
    
    res <- F_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Fsq,
                          f01,
                          f26)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7k) Bpa option
  #-------------------------------------------------------------------------------
  
  if("bpa" %in% stf.options){
    
    caseName <- "bpa"
    
    res <- B_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Bpa)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7l) Blim option
  #-------------------------------------------------------------------------------
  
  if("blim" %in% stf.options){
    
    caseName <- "blim"
    
    res <- B_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$Blim)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  #-------------------------------------------------------------------------------
  # 7m) MSYBtrigger option
  #-------------------------------------------------------------------------------
  
  if("MSYBtrigger" %in% stf.options){
    
    caseName <- "MSYBtrigger"
    
    res <- B_scaling_fun( stf,
                          FuY,
                          CATCH,
                          RECS,
                          referencePoints,
                          referencePoints$MSYBtrigger)
    
    stf.table[caseName,"Fbar 2-6 A",]                                 <- quantMeans(res$stf@harvest[f26,FcY,"A"])
    stf.table[caseName,grep("Fbar 0-1 ",dimnames(stf.table)$values),] <- aperm(quantMeans(res$stf@harvest[f01,FcY,c("B","D")]),c(2:6,1))
    stf.table[caseName,"Fbar 1-3 C",]                                  <- aperm(quantMeans(res$stf@harvest[f13,FcY,c("C")]),c(2:6,1))
    stf.table[caseName,"Fbar 2-6",]                                   <- quantMeans(unitSums(res$stf@harvest[f26,FcY,]))
    stf.table[caseName,"Fbar 0-1",]                                   <- quantMeans(unitSums(res$stf@harvest[f01,FcY,]))
    stf.table[caseName,grep("Catch ",dimnames(stf.table)$values),]    <- res$stf@catch[,FcY]
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[1],]    <- res$ssb.FcY
    stf.table[caseName,grep("SSB",dimnames(stf.table)$values)[2],]    <- res$ssb.CtY
  
    stf.res <- bind_rows(stf.res,
                         tibble(scenario = caseName, stf=list(res$stf)))
  }
  
  ## -----------------------------------------------------------------------------
  ## 3  medium term projections (NEED TO CHECK THE FMSYAR RULE IMPLEMENTATION: SEEMS WRONG HERE)
  ## -----------------------------------------------------------------------------
  
  # First run the Fmsy AR again
  res2 <- fmsyAR_fun_notransfer(stf,
                                FuY,
                                TACS,
                                RECS,
                                referencePoints,
                                TAC_var,
                                f01,
                                f26)
  
  # or the Fsq version
  # res2 <- F_scaling_fun( stf,
  #                       FuY,
  #                       CATCH,
  #                       RECS,
  #                       referencePoints,
  #                       referencePoints$Fsq,
  #                       f01,
  #                       f26)
  
  # or the F=0.2 version
  # res2 <- F_scaling_fun( stf,
  #                        FuY,
  #                        CATCH,
  #                        RECS,
  #                        referencePoints,
  #                        0.2,
  #                        f01,
  #                        f26)
  
  # set up the stf2 object  
  stf2  <- window(res2[["stf"]],start=an(dms$year)[1],end=(4+rev(an(dms$year))[1]))
  FuY2  <- ac((an(CtY)+1):(an(CtY)+4))
  stf2@harvest[,CtY]  <- stf2@harvest[,FcY]
  stf2@harvest[,FuY2] <- stf2@harvest[,FcY]
  
  # fill the slots
  for(i in c("stock.wt", "catch.wt", "mat", "m", "harvest.spwn", "m.spwn")){
    slot(stf2,i)[,FuY2] <- slot(stf2,i)[,CtY]
  }
  
  # fill the recruitment
  stf2@stock.n["0",FuY2] <- stf2@stock.n["0",CtY]
  
  # calculate catch in CtY
  for(i in dms$unit){
    stf2@catch.n[,CtY,i]     <- stf2@stock.n[,CtY,i]*(1-exp(-unitSums(stf2@harvest[,CtY])-stf2@m[,CtY,i]))*(stf2@harvest[,CtY,i]/(unitSums(stf2@harvest[,CtY])+stf2@m[,CtY,i]))
    stf2@catch[,CtY,i]       <- computeCatch(stf2[,CtY,i])
  }
  
  # calculate stock and catch in future years
  for(y in FuY2) {
    j <- an(y)
    
    # age 2 to plusgroup-1
    stf2@stock.n[2:(dims(stf2)$age-1),y]   <- 
      (stf2@stock.n[,ac(j-1),1]*exp(-unitSums(stf2@harvest[,ac(j-1),])-stf2@m[,ac(j-1),1]))[ac(range(stf2)["min"]:(range(stf2)["max"]-2))]
    
    # plusgroup
    stf2@stock.n[dims(stf2)$age,y]         <- 
      apply((stf2@stock.n[,ac(j-1),1]*
               exp(-unitSums(stf2@harvest[,ac(j-1)])-stf2@m[,ac(j-1),1]))[ac((range(stf2)["max"]-1):range(stf2)["max"]),],2:6,sum,na.rm=T)
    #ssb
    stf2@stock[,y] <- quantSums( stf2@stock.n[,y,1] * stf2@stock.wt[,y,1] *
                                   exp(-unitSums(stf2@harvest[,y])*stf2@harvest.spwn[,y,1]-stf2@m[,y,1] *
                                         stf2@m.spwn[,y,1]) * stf2@mat[,y,1])[,y,1]
    
    
    for(i in dms$unit){
      stf2@catch.n[,y,i]     <- stf2@stock.n[,y,i]*(1-exp(-unitSums(stf2@harvest[,y])-stf2@m[,y,i]))*(stf2@harvest[,y,i]/(unitSums(stf2@harvest[,y])+stf2@m[,y,i]))
      stf2@catch[,y,i]       <- computeCatch(stf2[,y,i])
    }
  }
  
  # calculate stock size
  stf2@stock[,] <- quantSums( stf2@stock.n[,,1] * stf2@stock.wt[,,1] *
                                exp(-unitSums(stf2@harvest[,])*stf2@harvest.spwn[,,1]-stf2@m[,,1] *
                                      stf2@m.spwn[,,1]) * stf2@mat[,,1])[,,1]
  
  
  
  # convert to dataframe
  stf2.df <- as.data.frame(stf2) %>% 
    #   bind_rows(., mutate(as.data.frame(fbar(stf2)[,,1]), slot="FA2-6")) %>% 
    bind_rows(., mutate(as.data.frame(quantMeans(unitSums(stf2@harvest[f26,,1]))), slot="FA2-6")) %>% 
    bind_rows(., mutate(as.data.frame(rec(stf2)[,,1]), slot="rec", age=as.character(age)))
  
  # stf2.df %>% 
  #   filter(slot=="FA2-6") %>% 
  #   print()
  
  ## -----------------------------------------------------------------------------
  ## 4  saving outputs
  ## -----------------------------------------------------------------------------
  
  if(TAC_CD_advice){
    save(stf.res, stf.table, res2, stf2, stf2.df, 
         file=file.path('.','model',paste0('forecast_',assessment_year,'.RData')))
  }else{
    save(stf.res, stf.table, res2, stf2, stf2.df, 
         file=file.path('.','model',paste0('forecast_',assessment_year,'_CDsq.RData')))
  }
}


