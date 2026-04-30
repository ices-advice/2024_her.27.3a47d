## Prepare plots and tables for report

## Before: forecast.csv (output)
## After:  forecast.csv (report)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(doBy)
library(ggplot2)
library(dplyr)
library(ggplotFL)   # install.packages("ggplotFL", repos="http://flr-project.org/R")
library(tidyverse)
library(reshape2)

mkdir(file.path('./',"report",'nego'))

data.source <- file.path('.','bootstrap','data')
model.dir  <- file.path('.','model')
report.dir  <- file.path('.','report','nego')

sourceDir("bootstrap/initial/software/utilities")


## -----------------------------------------------------------------------------
## plotting individual stfs
## -----------------------------------------------------------------------------
assessment_year     <- '2024'

# load forecast results (stf & stf2.df)
load(file.path(model.dir,paste0('forecast_',assessment_year,'_nego.RData')))

load(file.path('.','bootstrap','data',
               paste0('NSAS_HAWG',assessment_year,'_sf.RData')))

# load TAC
TACTab            <- read.taf("bootstrap/data/stf_NSAS_TAC.csv")
adviceTab         <- read.taf("bootstrap/data/stf_NSAS_advice.csv")
transferTab         <- read.taf("bootstrap/data/stf_NSAS_transfer.csv")
catchTab         <- read.taf("bootstrap/data/stf_NSAS_catches.csv")
catchTab$WBSS <- catchTab$C_WBSS+catchTab$D_WBSS+catchTab$`22-24_WBSS`+catchTab$IVaE_WBSS
catchTab$NSAS <- catchTab$A+catchTab$B+catchTab$C_NSAS+catchTab$D_NSAS
catchTab <- left_join(catchTab,transferTab,by='year')
catchTab <- catchTab %>% select(c('year','NSAS','WBSS','C.prop','D.prop'))

adviceTab <- pivot_longer(adviceTab,!year,names_to='fleets',values_to='data')
adviceTab$type <- 'advise'
TACTab <- pivot_longer(TACTab,!year,names_to='fleets',values_to='data')
TACTab$type <- 'TAC'
adviseTAC <- rbind(adviceTab,TACTab)

catchAssessment <- as.data.frame(NSH@catch)

NSH.df <- 
  as.data.frame(NSH) %>% 
  bind_rows(as.data.frame(ssb(NSH)) %>% mutate(slot="ssb")) %>% 
  bind_rows(as.data.frame(rec(NSH)) %>% mutate(slot="rec", age=as.character(age))) %>% 
  bind_rows(as.data.frame(fbar(NSH)) %>% mutate(slot="fbar")) %>% 
  mutate(source="assessment")

##################################################
# plotting
##################################################

adviseTAC.tot <- subset(adviseTAC,fleets != 'F') %>% group_by(type,year) %>% summarize(data=sum(data,na.rm = T))

taf.png(file.path(report.dir,"catch_TAC"))
print(ggplot(subset(adviseTAC.tot,year >= 2018 & year != 2025),
       aes(x=as.factor(year),y=data,fill=type))+
  theme_bw()+
  ylab('Catch')+
  geom_bar(stat = 'identity',position="dodge"))
dev.off()


