altRunsFlag   <- F
refPointsFlag <- F
profilingFlag <- F

source("model_sf.R")
source("model_mf.R")

if(altRunsFlag == T){
  #################################
  # running aleternative sf runs
  #################################
  source('model_sf_altRuns.R')
  
  #################################
  # running aleternative mf runs
  #################################
  source('model_mf_altRuns.R')
}

if(refPointsFlag == T){
  #################################
  # ref points as of HAWG2024
  #################################
  rm(list=ls())
  filename  = "NSAS_HAWG2024_sf.RData"
  source('model_refpoints.r')
  
  #################################
  # ref points as of IBP2021
  #################################
  rm(list=ls())
  filename  = "NSAS_HAWG2024_sf_IBP2021.RData"
  source('model_refpoints.r')
  
  ################################################
  # retrospective ref points with HAWG2024 config
  ################################################
  rm(list=ls())
  filename  = "NSAS_HAWG2024_sf_retro.RData"
  source('model_refpoints_retro.r')
  
  ################################################
  # retrospective ref points with IBP config
  ################################################
  rm(list=ls())
  filename  = "NSAS_HAWG2024_sf_SMS2019_2_retro.RData"
  source('model_refpoints_retro.r')
}


if(profilingFlag == T){
  ################################################
  # assessment profiling
  ################################################
  source('model_scanM.R')
  source('model_scanM_peels_sf.R')
  source('model_scanM_SMS2019.R')
  source('model_scanM_smsAll_sf.R')
}
