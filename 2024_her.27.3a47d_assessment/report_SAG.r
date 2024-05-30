# ============================================================================
# UPload Assessment Data into SAG
#
# 17/08/2021 Updated for HAWG 2021
# 30/04/2022 Updated for HAWG 2022
# ============================================================================

rm(list=ls())

library(icesTAF)

taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)

library(tidyverse)
library(icesSAG)  # install.packages("icesSAG")

# Warning in install.packages :
#   package ‘icesSAG’ is not available for this version of R
# 
# A version of this package for your version of R might be available elsewhere,
# see the ideas at
# https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

runName <- 'NSAS_HAWG2024'

model.save      <- file.path(".",'model')
load(file.path(model.save,paste0(runName,'_sf.RData')))

# use token
options(icesSAG.use_token = TRUE)


# Set years and ranges
FiY   <- dims(NSH)$minyear
DtY   <- dims(NSH)$maxyear
LaY   <- dims(NSH)$maxyear+1
nyrs  <- ((DtY)-(FiY))+1
nyrs2 <- ((LaY)-(FiY))+1

# Meta information
stockkeylabel  <- "her.27.3a47d"
assessmentyear <- 2024
contactperson  <- "benoit.berges@wur.nl"

# SSB in intermediate year 
SSBint <- 1386241   # to be updated after the forecast

# Create the input data for uploading  
info     <- icesSAG::stockInfo(StockCategory = 1,
                               ModelType='A',
                               Purpose = "Advice",
                               ModelName = 'SAM',
                               StockCode      = stockkeylabel, 
                               AssessmentYear = assessmentyear, 
                               ContactPerson  = contactperson)

info$StockCategory             <- "1"
info$MSYBtrigger               <- 1130070
info$Blim                      <- 828874
info$Bpa                       <- 903707
info$Flim                      <- 0.39
info$Fpa                       <- 0.32
info$FMSY                      <- 0.32
info$Fage                      <- "2-6" 
info$RecruitmentAge            <- 0
info$CatchesLandingsUnits      <- "t"
info$RecruitmentDescription    <- "wr"
info$RecruitmentUnits          <- "NE3" 
info$FishingPressureDescription<- "F"
info$FishingPressureUnits      <- NA 
info$StockSizeDescription      <- "SSB"
info$StockSizeUnits            <- "t"
info$Purpose                   <- "Advice"
info$CustomSeriesName1         <- "model catch"
info$CustomSeriesName2         <- "model catch low"
info$CustomSeriesName3         <- "model catch high"
info$CustomSeriesName4         <- "F0-1"
info$CustomSeriesName5         <- "F2-6"
info$CustomSeriesName6         <- "F7-8"
info$CustomSeriesUnits1        <- "t"
info$CustomSeriesUnits2        <- "t"
info$CustomSeriesUnits3        <- "t"
info$CustomSeriesUnits4        <- NA
info$CustomSeriesUnits5        <- NA
info$CustomSeriesUnits6        <- NA
info$ModelName                 <- "SAM"
info$ModelType                 <- "A"
info$ConfidenceIntervalDefinition <- "95%"

# Create the fish data
fishdata                          <- stockFishdata(FiY:LaY)

fishdata$Catches[1:nyrs]          <- an(NSH@landings)[1:nyrs]

fishdata$Low_Recruitment          <- rec(NSH.sam)$lbnd
fishdata$Recruitment              <- rec(NSH.sam)$value
fishdata$High_Recruitment         <- rec(NSH.sam)$ubnd 

fishdata$Low_StockSize[1:nyrs]    <- ssb(NSH.sam)$lbnd[1:nyrs]
fishdata$StockSize                <- c(ssb(NSH.sam)$value[1:nyrs], SSBint)
fishdata$High_StockSize[1:nyrs]   <- ssb(NSH.sam)$ubnd[1:nyrs]

fishdata$Low_TBiomass[1:nyrs]     <- tsb(NSH.sam)$lbnd[1:nyrs]
fishdata$TBiomass                 <- tsb(NSH.sam)$value
fishdata$High_TBiomass[1:nyrs]    <- tsb(NSH.sam)$ubnd[1:nyrs]

fishdata$Low_FishingPressure[1:nyrs] <- fbar(NSH.sam)$lbnd[1:nyrs]
fishdata$FishingPressure[1:nyrs]     <- fbar(NSH.sam)$value[1:nyrs]
fishdata$High_FishingPressure[1:nyrs]<- fbar(NSH.sam)$ubnd[1:nyrs]

fishdata$CustomSeries1[1:nyrs]    <- catch(NSH.sam)$value[1:nyrs]
fishdata$CustomSeries2[1:nyrs]    <- catch(NSH.sam)$lbnd[1:nyrs]
fishdata$CustomSeries3[1:nyrs]    <- catch(NSH.sam)$ubnd[1:nyrs]

fishdata$CustomSeries4[1:nyrs]    <- c(quantMeans(harvest(NSH.sam)[ac(0:1),]))[1:nyrs]
fishdata$CustomSeries5[1:nyrs]    <- c(quantMeans(harvest(NSH.sam)[ac(2:6),]))[1:nyrs]
fishdata$CustomSeries6[1:nyrs]    <- c(quantMeans(harvest(NSH.sam)[ac(7:8),]))[1:nyrs]

# View(fishdata)

# upload to SAG
key <- icesSAG::uploadStock(info, fishdata)

# Get SAG settings
# getSAGSettingsForAStock(assessmentKey=key) %>% View()
# getSAGSettingsForAStock(assessmentKey=14095) %>% View()

# Add comment to SAG settings
# setSAGSettingForAStock(assessmentKey=key, 
#                        chartKey=0,
#                        settingKey=21,
#                        settingValue="My text for the comment field",
#                        copyNextYear=FALSE) 

# plot F's
# fishdata %>%
#   dplyr::select(Year, FishingPressure, CustomSeries4, CustomSeries6) %>%
#   setNames(c("year", "F26", "F01", "F78")) %>%
#   gather(key=F, value=value, F26:F78) %>%
#   filter(year >= 1980) %>%
# 
#   ggplot(aes(x=year, y=value, group=F)) +
#   theme_bw() +
#   geom_line(aes(colour=F), size=1)

