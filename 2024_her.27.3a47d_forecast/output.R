## Extract results of interest, write TAF output tables

## Before: forecast.RData, transcript.txt (model)
## After:  forecast.csv, forecast.RData, transcript.txt (output)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)

mkdir("output")

## Copy forecast results to output folder
#cp("model/*", "output")

assessment_year     <- '2024'

data.source <- file.path('.','bootstrap','data')

for(TAC_CD_advice in c(TRUE,FALSE)){
  ## Write forecast table to file
  
  if(TAC_CD_advice){
    load(file.path('.','model',paste0('forecast_',assessment_year,'.RData')))
  }else{
    load(file.path('.','model',paste0('forecast_',assessment_year,'_CDsq.RData')))
  }
  
  forecast <- xtab2taf(stf.table[,,"50%"])
  names(forecast) <- c("Basis", "Fbar26A", "Fbar01B", "Fbar13C", "Fbar01D",
                       "Fbar26", "Fbar01", "CatchA", "CatchB", "CatchC", "CatchD",
                       "SSB1", "SSB2")
  
  if(TAC_CD_advice){
    write.taf(forecast, dir="output",file = paste0('forecast_',assessment_year,'.csv'))
  }else{
    write.taf(forecast, dir="output",file = paste0('forecast_',assessment_year,'_CDsq.csv'))
  }
}