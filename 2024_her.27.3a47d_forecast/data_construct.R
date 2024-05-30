## Preprocess data, write TAF data tables

# Code to do the multifleet short term forecast for North Sea Herring

rm(list=ls())

library(icesTAF)

taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)   # rec

data.source <- file.path('.','bootstrap','data')

assessment_year   <- '2024'

### ============================================================================
### Construct packages object
### ============================================================================

library(minpack.lm)  # nls.lm, in fleet.harvest
sourceDir("bootstrap/initial/software/utilities")

mkdir("data")

referencePoints <- list(Fmsy = 0.32,
                        Fsq  = NA,
                        Flim = 0.39,
                        Fpa  = 0.33,
                        Blim = 828874,
                        Bpa  = 903707,
                        MSYBtrigger = 1131601,
                        Ftarget  = 0, # set this for individual cases
                        F01   = 0.05,
                        Btrigger = 0)

# flag on TAC assumptions for C and D fleet. 
# If true, one takes TAC from WBSS advice
# If false, sq in TAC, i.e. one takes TAC from ImY (fed from WBSS advice)
#TAC_CD_advice   <- TRUE

for(TAC_CD_advice in c(TRUE,FALSE)){
  
  ## -----------------------------------------------------------------------------
  ## 1  Load data and set location to data folder
  ## -----------------------------------------------------------------------------
  
  load(file.path('.','bootstrap','data',
                 paste0('NSAS_HAWG',assessment_year,'_sf.RData')))
  load(file.path('.','bootstrap','data',
                 paste0('NSAS_HAWG',assessment_year,'_mf.RData')))
  
  #-------------------------------------------------------------------------------
  # 2) setup control variables
  #-------------------------------------------------------------------------------
  
  DtY   <- ac(range(NSH)["maxyear"]) #Data year
  ImY   <- ac(an(DtY)+1)             #Intermediate year
  FcY   <- ac(an(DtY)+2)             #Forecast year
  CtY   <- ac(an(DtY)+3)             #Continuation year
  CtY1  <- ac(an(DtY)+4)
  FuY   <- c(ImY,FcY,CtY)            #Future years
  
  # slots copied from last year of assessment
  yrs1        <- list("m.spwn","harvest.spwn","stock.wt")
  
  # slots copied as mean of last 3 years of assessment
  yrs3        <- list("mat")
  
  # slots copied as mean of last 5 years of assessment
  yrs5        <- list("m")
  
  dsc         <- "North Sea Herring"
  nam         <- "NSH"
  dms         <- dimnames(NSH@m)
  dms$year    <- c(rev(rev(dms$year)[1:3]),ImY,FcY,CtY)
  dms$unit    <- c("A","B","C","D")
  
  f01         <- ac(0:1)
  f13         <- ac(1:3)
  f26         <- ac(2:6)
  
  stf.options <- c('fmsyAR',
                   'fmsyAR_Btarget',
                   'fmsyAR_TACrule_transfer',
                   'fmsyAR_TACrule_notransfer',
                   "fmsyAR_notransfer",
                   "fmsyAR_Btarget_notransfer",
                   "fmsy",
                   "nf",
                   "tacro",
                   "-15%",
                   "+15%",
                   "fsq",
                   "fpa",
                   "flim",
                   "bpa",
                   "blim",
                   "MSYBtrigger")
  # mp    =according to management plan,
  # +/-%  =TAC change,
  # nf    =no fishing,
  # bpa   =reach BPA in CtY or if ssb > bpa fish at fpa,
  # tacro =same catch as last year
  # fmsy  = fmsy implemented
  # flim,fpa,blim,bpa = reach F or biomass targets
  
  #-------------------------------------------------------------------------------
  # 3) TAC information, reference points and management rule
  #
  # Note 1: TACS are the realised NSAS catches for the different fleets
  # Note 2: TACS.orig are the original set TACs
  #-------------------------------------------------------------------------------
  
  ############## TAC variables in the ImY ##############
  # initialize TAC variables array
  TAC_var                   <- list(CtransferFcY = NA,
                                    CtransferImY = NA,
                                    DtransferImY = NA,
                                    Csplit = NA,
                                    Dsplit = NA,
                                    WBSS_NSAS = NA,
                                    Buptake = NA,
                                    Duptake = NA)
  
  # reading catches table
  catchAreaTab            <- read.taf("bootstrap/data/stf_NSAS_catches.csv") # historical C fleet WBSS/NSAS split (Henrik Mosegaard)
  rownames(catchAreaTab) <- catchAreaTab$year
  
  # reading transfer variables table
  # Ctransfer. Transfer of TAC from IIIa to IVa for C fleet in assessment year
  # Note: C transfer in table is the percentage transferred to NSAS from WBSS
  transferTab            <- read.taf("bootstrap/data/stf_NSAS_transfer.csv")
  rownames(transferTab) <- transferTab$year
  transferTab$D_TAC.NS <- NA
  
  # TAC information for current year
  TACTab            <- read.taf("bootstrap/data/stf_NSAS_TAC.csv")
  rownames(TACTab)  <- TACTab[,1]
  
  # tAC advice information
  adviceTab             <- read.taf("bootstrap/data/stf_NSAS_advice.csv")
  rownames(adviceTab)   <- adviceTab[,1]
  
  # proportion of WBSS in the C and D fleets
  catchAreaTab$C.split    <- catchAreaTab$C_WBSS/(catchAreaTab$C_NSAS+catchAreaTab$C_WBSS)
  catchAreaTab$D.split    <- catchAreaTab$D_WBSS/(catchAreaTab$D_NSAS+catchAreaTab$D_WBSS)
  catchAreaTab$NS.split   <- catchAreaTab$IVaE_WBSS/(catchAreaTab$A+catchAreaTab$IVaE_WBSS)
  catchAreaTab$totalEU3a  <- catchAreaTab$C_WBSS+catchAreaTab$C_NSAS-catchAreaTab$C_NO+catchAreaTab$D_NSAS+catchAreaTab$D_WBSS
  catchAreaTab$total3a    <- catchAreaTab$C_WBSS+catchAreaTab$C_NSAS+catchAreaTab$D_NSAS+catchAreaTab$D_WBSS
  catchAreaTab$CshareEU3a <- (catchAreaTab$C_WBSS+catchAreaTab$C_NSAS-catchAreaTab$C_NO)/catchAreaTab$totalEU3a
  catchAreaTab$DshareEU3a <- 1-catchAreaTab$CshareEU3a
  catchAreaTab$Cuptake    <- (catchAreaTab$C_WBSS+catchAreaTab$C_NSAS)/TACTab[ac(catchAreaTab$year),]$C
  catchAreaTab$Duptake    <- (catchAreaTab$D_WBSS+catchAreaTab$D_NSAS)/TACTab[ac(catchAreaTab$year),]$D
  catchAreaTab$Buptake    <- catchAreaTab$B/TACTab[ac(catchAreaTab$year),]$B
  catchAreaTab$BDuptake   <- catchAreaTab$B/(TACTab[ac(catchAreaTab$year),]$B+TACTab[ac(catchAreaTab$year),]$D*transferTab[DtY,]$D.prop)
  transferTab[ImY,]$C_TAC.NS  <- TACTab[ImY,]$C-(transferTab[ImY,]$CD_TAC_EU*catchAreaTab[DtY,]$CshareEU3a+catchAreaTab[ImY,]$C_NO)
  transferTab[ImY,]$D_TAC.NS  <- transferTab[ImY,]$CD_TAC_EU*catchAreaTab[DtY,]$DshareEU3a

  # C split as average over the last 3 years
  TAC_var$Csplit      <- mean(1-catchAreaTab[ac((an(DtY)-2):an(DtY)),]$C.split)    # Proportion NSAS in C fleet catch; 3 year average 
  TAC_var$Dsplit      <- mean(1-catchAreaTab[ac((an(DtY)-2):an(DtY)),]$D.split)    # Proportion NSAS in C fleet catch; 3 year average 
  TAC_var$WBSS_NSAS   <- mean(catchAreaTab[ac((an(DtY)-2):an(DtY)),]$NS.split)
  TAC_var$CshareEU3a  <- catchAreaTab[DtY,]$CshareEU3a
  TAC_var$DshareEU3a  <- catchAreaTab[DtY,]$DshareEU3a
  TAC_var$Buptake         <- catchAreaTab[DtY,]$BDuptake #mean(uptakeTab[ac((an(DtY)-2):an(DtY)),'B'])   # Uptake of Bfleet as 3 year average
  TAC_var$Duptake         <- catchAreaTab[DtY,]$Duptake #mean(catchAreaTab[ac((an(DtY)-2):an(DtY)),]$Duptake)   # Uptake of the Dfleet as 3 year average
  TAC_var$CtransferImY <- transferTab[ImY,]$C_TAC.NS
  TAC_var$CtransferFcY <- transferTab[ImY,]$C_TAC.NS/TACTab[ImY,]$C
  TAC_var$DtransferImY <- TACTab[ImY,]$D*transferTab[ImY,]$D.prop-transferTab[ImY,]$D_TAC.NS
  TAC_var$DtransferFcY <- transferTab[ImY,]$D_TAC.NS/TACTab[ImY,]$D
  
  ############## TAC for the different fleets in the intermediate year ##############
  # TAC A is TACNSA; HER/4AB. + HER/4CXB7D 
  # TAC B is TACNSB; HER/2A47DX 
  # TAC C is TAC3aC; HER/03A.
  # TAC D is TAC3aD; HER/03A-BC
  
  # Create TAC objects for current year, forecast year and last year
  TACS        <- FLQuant(NA,dimnames=list(age="all",year=FuY,unit=c("A","B","C","D","F"),
                                          season="all",area=1,iter=1:dims(NSH)$iter))
  
  TACS[,ImY,'A'] <- TACTab[ImY,'A']
  TACS[,ImY,'B'] <- TACTab[ImY,'B']
  # assume 0 TAC for FcY and CtY for C and D fleet
  TACS[,ImY,'C'] <- TACTab[ImY,'C']
  TACS[,ImY,'D'] <- TACTab[ImY,'D']
  
  # TAC C and D fleet in FcY and CtY
  # set TAC to 0.1 if 0 catch (for optimizers)
  if(TAC_CD_advice == TRUE){
    TACS[,c(FcY,CtY),'C'] <- adviceTab[FcY,'C']
    if(TACS[,FcY,'C'] == 0){
      TACS[,c(FcY,CtY),'C'] <- 0.1
    }
  }else{
    TACS[,c(FcY,CtY),'C'] <- TACTab[ImY,'C']
  }
  
  # D fleet
  if(TAC_CD_advice == TRUE){
    TACS[,c(FcY,CtY),'D'] <- adviceTab[FcY,'D']
    if(TACS[,FcY,'D'] == 0){
      TACS[,c(FcY,CtY),'D'] <- 0.1
    }
  }else{
    TACS[,c(FcY,CtY),'D'] <- TACTab[ImY,'D']
  }
  
  TACS[,c(ImY,FcY,CtY),'F'] <- TACTab[c(ImY,FcY,FcY),'F']
  
  ############## realised catches ##############
  CATCH        <- FLQuant(NA,dimnames=list( age="all",year=FuY,unit=c("A","B","C","D"),
                                            season="all",area=1,iter=1:dims(NSH)$iter))
  
  CATCH[,ImY,'A'] <- TACS[,ImY,'A'] + TAC_var$CtransferImY-(TACS[,ImY,'A'] + TAC_var$CtransferImY)*TAC_var$WBSS_NSAS
  #CATCH[,ImY,'A'] <- TACS[,ImY,'A'] + TAC_var$Ctransfer*TACS[,ImY,'C'] - (TACS[,ImY,'A'] + TAC_var$Ctransfer*TACS[,ImY,'C'])*TAC_var$WBSS_NSAS
  CATCH[,ImY,'B'] <- (TACS[,ImY,'B']+TAC_var$DtransferImY)*TAC_var$Buptake
  CATCH[,ImY,'C'] <- (transferTab[ImY,]$CD_TAC_EU*catchAreaTab[DtY,]$CshareEU3a+catchAreaTab[ImY,]$C_NO)*TAC_var$Csplit
  CATCH[,ImY,'D'] <- transferTab[ImY,]$CD_TAC_EU*TAC_var$DshareEU3a*TAC_var$Dsplit
  
  # assume same catch as ImY
  # This is an alternative to the previous procedure of estimating the B-fleet catch from the Management plan
  #CATCH[,c(FcY,CtY),'B'] <- (TACS[,ImY,'B']+TACS[,FcY,'D']*TAC_var$DtransferImY)*TAC_var$Buptake
  
  # zero catch in FcY and CtY for C and D fleet, because of zero catch advice for WBSS herring
  if(TACS[,FcY,'C'] == 0.1){
    CATCH[,c(FcY,CtY),'C'] <- 0.1
  }else{
    CATCH[,c(FcY,CtY),'C'] <- TACS[,FcY,'C']*(1-TAC_var$CtransferFcY)*TAC_var$Csplit
  }
  
  if(TACS[,FcY,'D'] == 0.1){
    CATCH[,c(FcY,CtY),'D'] <- 0.1
  }else{
    # assume no fishing
    CATCH[,c(FcY,CtY),'D'] <- TACS[,FcY,'D']*(1-TAC_var$DtransferImY)*TAC_var$Dsplit*TAC_var$Duptake
  }
  
  #-------------------------------------------------------------------------------
  # 4) Recruitment for intermediate, forecast and continuation years
  #    ImY recruitment: comes from sf assessment (informed by IBTS0)
  #    FcY and CtY: weighted average over 10 year using uncertainty as weighting
  #-------------------------------------------------------------------------------
  # Retrieve uncertainty estimate on recruitment estimates by the model for weighted average
  recWeights  <- subset(NSH.sam@params, name=="logR")$std.dev^2
  
  # fill in recruitment for the different years
  RECS <- FLQuants( "ImY" =FLQuant(subset(rec(NSH.sam),year==ImY)$value),
                    "FcY" =exp(apply(log(rec(NSH)[,ac((an(DtY)-9):DtY)]),3:6,weighted.mean,w=1/rev(rev(recWeights)[2:11]),na.rm=T)),
                    "CtY" =exp(apply(log(rec(NSH)[,ac((an(DtY)-9):DtY)]),3:6,weighted.mean,w=1/rev(rev(recWeights)[2:11]),na.rm=T)))
  
  #RECS <- FLQuants( "ImY" =FLQuant(subset(rec(NSH.sam),year==ImY)$value),
  #                  "FcY" =FLQuant(subset(rec(NSH.sam),year==ImY)$value),
  #                  "CtY" =FLQuant(subset(rec(NSH.sam),year==ImY)$value))
  
  #-------------------------------------------------------------------------------
  # 5) Setup stock file
  #-------------------------------------------------------------------------------
  
  # Create an stf object 
  stf         <- FLStock(name=nam,
                         desc=dsc,
                         m=FLQuant(NA,
                                   dimnames=dms))
  
  range(stf)['maxfbar'] <- 6
  range(stf)['minfbar'] <- 2
  
  # copy fields from current assessment to all fleets
  for(i in dms$unit)
    stf[,,i]  <- window(NSH,start=an(dms$year)[1],end=rev(an(dms$year))[1])
  
  units(stf)  <- units(NSH)
  
  # Fill slots that are the same for all fleets
  for(i in c(unlist(yrs1),unlist(yrs3),unlist(yrs5))){
    if(i %in% unlist(yrs1)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- slot(NSH,i)[,DtY]}}
    if(i %in% unlist(yrs3)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(NSH,i)[,ac((an(DtY)-2):an(DtY))])}}
    if(i %in% unlist(yrs5)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(NSH,i)[,ac((an(DtY)-4):an(DtY))])}}
  }
  
  # taking data from multifleet assessment
  for(i in 1:length(dms$unit)){
    
    # Fill the D fleet with the data from the B fleet; same selection pattern
    if(i == 4) {
      # i = 4
      stf@harvest[,ac(an(dms$year)[1]:(an(dms$year)[1]+2)),4]     <- NSH3f.sam@harvest[,ac(an(dms$year)[1]:(an(dms$year)[1]+2)),,,2]
      stf@harvest[,FuY,4]     <- NSH3f.sam@harvest[,ImY,,,2]
      stf@catch.wt[,FuY,4]    <- yearMeans(NSH3f@catch.wt[,dms$year[1:3],,,2])
      stf@landings.wt[,FuY,4] <- yearMeans(NSH3f@catch.wt[,dms$year[1:3],,,2])
    }
    
    if(i != 4) {
      stf@harvest[,ac(an(dms$year)[1]:(an(dms$year)[1]+2)),i]     <- NSH3f.sam@harvest[,ac(an(dms$year)[1]:(an(dms$year)[1]+2)),,,i]
      stf@harvest[,FuY,i]     <- NSH3f.sam@harvest[,ImY,,,i]
      stf@catch.wt[,FuY,i]    <- yearMeans(NSH3f@catch.wt[,dms$year[1:3],,,i])
      stf@landings.wt[,FuY,i] <- yearMeans(NSH3f@catch.wt[,dms$year[1:3],,,i])
    }
  }

  # Fill slots that have no meaning for NSAS
  stf@discards.n[]          <- 0
  stf@discards[]            <- 0
  stf@discards.wt[]         <- 0
  
  #-------------------------------------------------------------------------------
  # 5) Compute F in intermediate year
  #-------------------------------------------------------------------------------
  
  # stock number in ImY, using SAM estimate
  for(i in dms$unit)
    stf@stock.n[,ImY,i]     <- NSH.sam@stock.n[,ImY]
  
  # computing F 
  stf@harvest[,ImY]         <- fleet.harvest(stk=stf,
                                             iYr=ImY,
                                             CATCH=CATCH[,ImY])
  for(i in dms$unit){
    stf@catch.n[,ImY,i]     <- stf@stock.n[,ImY,i]*(1-exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,i]))*(stf@harvest[,ImY,i]/(unitSums(stf@harvest[,ImY])+stf@m[,ImY,i]))
    stf@catch[,ImY,i]       <- computeCatch(stf[,ImY,i])
    stf@landings.n[,ImY,i]  <- stf@catch.n[,ImY,i]
    stf@landings[,ImY,i]    <- computeLandings(stf[,ImY,i])
  }
  
  
  #Intermediate year stf option table
  stf.table <- array(NA,dim=c(c(length(stf.options)+1),12,3),
                     dimnames=list("options"=c("intermediate year",stf.options),
                                   "values" =c("Fbar 2-6 A","Fbar 0-1 B","Fbar 1-3 C","Fbar 0-1 D","Fbar 2-6","Fbar 0-1","Catch A","Catch B","Catch C","Catch D","SSB","SSB"),#,"C transfer"
                                   "stats"  =c("5%","50%","95%")))
  
  stf.table[1,"Fbar 2-6 A",]                                  <- quantMeans(stf@harvest[f26,ImY,"A"])
  stf.table[1,grep("Fbar 0-1 ",dimnames(stf.table)$values),]  <- aperm(quantMeans(stf@harvest[f01,ImY,c("B","D")]),c(2:6,1))
  stf.table[1,"Fbar 1-3 C",]                                  <- aperm(quantMeans(stf@harvest[f13,ImY,c("C")]),c(2:6,1))
  stf.table[1,"Fbar 2-6",]                                    <- quantMeans(unitSums(stf@harvest[f26,ImY,]))
  stf.table[1,"Fbar 0-1",]                                    <- quantMeans(unitSums(stf@harvest[f01,ImY,]))
  stf.table[1,grep("Catch ",dimnames(stf.table)$values),]     <- aperm(harvestCatch(stf,ImY),c(2:6,1))
  stf.table[1,grep("SSB",dimnames(stf.table)$values)[1],]     <- quantSums(stf@stock.n[,ImY,1] * stf@stock.wt[,ImY,1] *
                                                                             exp(-unitSums(stf@harvest[,ImY])*stf@harvest.spwn[,ImY,1]-stf@m[,ImY,1]*stf@m.spwn[,ImY,1]) *
                                                                             stf@mat[,ImY,1])
  
  referencePoints$Fsq <- stf.table["intermediate year","Fbar 2-6","50%"]
  
  
  #-------------------------------------------------------------------------------
  # 6) propagate stock number to forecast year setup
  #-------------------------------------------------------------------------------
  
  for(i in dms$unit) stf@stock.n[1,FcY,i]                     <- RECS$FcY
  for(i in dms$unit) stf@stock.n[2:(dims(stf)$age-1),FcY,i]   <- (stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac(range(stf)["min"]:(range(stf)["max"]-2)),]
  for(i in dms$unit) stf@stock.n[dims(stf)$age,FcY,i]         <- apply((stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac((range(stf)["max"]-1):range(stf)["max"]),],2:6,sum,na.rm=T)
  
  
  if(TAC_CD_advice){
    save(CATCH, 
         CtY, 
         f01, 
         f26,
         f13,
         FcY, 
         FuY, 
         RECS, 
         referencePoints, 
         stf, 
         stf.options,
         stf.table, 
         TAC_var, 
         TACS,
         file=file.path('.','data',paste0('input_',assessment_year,'.RData')))
  }else{
    save(CATCH, 
         CtY, 
         f01, 
         f26,
         f13,
         FcY, 
         FuY, 
         RECS, 
         referencePoints, 
         stf, 
         stf.options,
         stf.table, 
         TAC_var, 
         TACS,
         file=file.path('.','data',paste0('input_',assessment_year,'_CDsq.RData')))
  }
}

