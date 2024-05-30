## Prepare plots and tables for report

## Before: forecast.csv (output)
## After:  forecast.csv (report)

rm(list=ls())

library(icesTAF)
library(data.table)

mkdir("report")

data.source <- file.path('.','bootstrap','data')
model.dir  <- file.path('.','model')

assessment_year     <- '2024'

# load assessment, TAC and advice values
load(file.path(data.source,paste0('NSAS_HAWG',assessment_year,'_sf.RData')))
TACTab                  <- read.taf("bootstrap/data/stf_NSAS_TAC.csv")
adviceTab               <- read.taf("bootstrap/data/stf_NSAS_advice.csv")
advice_HAWG             <- sum(adviceTab[adviceTab$year == an(assessment_year),c(2:5)])#adviceTab$A[adviceTab$year == an(assessment_year)]+adviceTab$B[adviceTab$year == an(assessment_year)]
TACA_HAWG               <- TACTab$A[TACTab$year == an(assessment_year)]

adviceTab               <- pivot_longer(adviceTab,!year,names_to='fleets',values_to='data')
adviceTab$type          <- 'advise'
TACTab                  <- pivot_longer(TACTab,!year,names_to='fleets',values_to='data')
TACTab$type             <- 'TAC'
adviseTAC               <- rbind(adviceTab,TACTab)

# load forecast results
load(file.path('model',paste0('forecast_',assessment_year,'.RData')))

stf <- stf.res[stf.res$scenario=="fmsyAR_notransfer",2]$stf[[1]]

catchStf  <- as.data.frame(stf@catch[,,'A'])
ssbStf    <- as.data.frame(ssb(stf[,,'A']))

order_rows <- c('fmsyAR_notransfer',
                'fmsy',
                'nf',
                'tacro',
                'fsq',
                'fpa',
                'flim',
                'bpa',
                'blim',
                'MSYBtrigger',
                'fmsyAR_Btarget_notransfer')

order_columns <- c('Basis',
                   'Fbar26A',
                   'Fbar01B',
                   'Fbar13C',
                   'Fbar01D',
                   'Fbar26',
                   'Fbar01',
                   'CatchA',
                   'CatchB',
                   'CatchC',
                   'CatchD',
                   'total_catch',
                   'SSB1',
                   'SSB2',
                   'SSB_change',
                   'TAC_change',
                   'advice_change')

ImY_columns <- c('Fbar26A',
                 'Fbar01B',
                 'Fbar13C',
                 'Fbar01D',
                 'Fbar26',
                 'Fbar01',
                 'CatchA',
                 'CatchB',
                 'CatchC',
                 'CatchD',
                 'SSB1')

## forecast (round)
forecast  <- read.taf(file.path('output',paste0('forecast_',assessment_year,'.csv')))
forecast  <- rnd(forecast, "Fbar", 3, grep=TRUE)
forecast  <- rnd(forecast, "Catch|SSB", grep=TRUE)
initColNames    <- colnames(forecast)
forecast  <- cbind( forecast,
                        rep(0,20,1),
                        rep(0,20,1),
                        rep(0,20,1),
                        rep(0,20,1))
colnames(forecast) <- c(initColNames,'total_catch','SSB_change','TAC_change','advice_change')
forecast$total_catch <- forecast$CatchA+forecast$CatchB+forecast$CatchC+forecast$CatchD
forecast$SSB_change <- rnd((forecast$SSB1-forecast$SSB1[1])/forecast$SSB1[1]*100,digits = 1)
forecast$TAC_change <- rnd((forecast$CatchA-TACA_HAWG)/TACA_HAWG*100,digits = 1)
forecast$advice_change <- rnd((forecast$total_catch-advice_HAWG)/advice_HAWG*100,digits = 1)

forecast_CDsq <- read.taf(file.path('output',paste0('forecast_',assessment_year,'_CDsq.csv')))
forecast_CDsq <- rnd(forecast_CDsq, "Fbar", 3, grep=TRUE)
forecast_CDsq <- rnd(forecast_CDsq, "Catch|SSB", grep=TRUE)
initColNames  <- colnames(forecast_CDsq)
forecast_CDsq  <- cbind(  forecast_CDsq,
                          rep(0,20,1),
                          rep(0,20,1),
                          rep(0,20,1),
                          rep(0,20,1))
colnames(forecast_CDsq) <- c(initColNames,'total_catch','SSB_change','TAC_change','advice_change')
forecast_CDsq$total_catch <- forecast_CDsq$CatchA+forecast_CDsq$CatchB+forecast_CDsq$CatchC+forecast_CDsq$CatchD
forecast_CDsq$SSB_change <- rnd((forecast_CDsq$SSB1-forecast_CDsq$SSB1[1])/forecast_CDsq$SSB1[1]*100,digits = 1)
forecast_CDsq$TAC_change <- rnd((forecast_CDsq$CatchA-TACA_HAWG)/TACA_HAWG*100,digits = 1)
forecast_CDsq$advice_change <- rnd((forecast_CDsq$total_catch-advice_HAWG)/advice_HAWG*100,digits = 1)
forecast_CDsq$Basis[forecast_CDsq$Basis == 'fmsyAR_TACrule_transfer'] <- 'fmsyAR_TACrule_transfer'
forecast_CDsq$Basis[forecast_CDsq$Basis == 'fmsyAR_TACrule_notransfer'] <- 'fmsyAR_TACrule_notransfer'

# forecast table
forecast_final <- forecast
forecast_final <- forecast_final[forecast_final$Basis %in% order_rows,]
forecast_final[forecast_final$Basis == 'tacro',] <- forecast_CDsq[forecast_CDsq$Basis == 'tacro',]
forecast_final <- forecast_final[match(order_rows,forecast_final$Basis),]
forecast_final <- rbind(forecast_final,forecast_CDsq[forecast_CDsq$Basis == 'fmsyAR_TACrule_transfer',])
forecast_final <- rbind(forecast_final,forecast_CDsq[forecast_CDsq$Basis == 'fmsyAR_TACrule_notransfer',])
forecast_final <- forecast_final[,match(order_columns,colnames(forecast_final))]

# ImY assumptions
ImY_assumptions <- data.table(  F26 = rnd(drop(unitSums(fbar(stf)[,assessment_year])),digits = 3),
                                SSBImY = rnd(drop(quantSums(stf@stock.n[,assessment_year,1] * stf@stock.wt[,assessment_year,1] *
                                                     exp(-unitSums(stf@harvest[,assessment_year])*stf@harvest.spwn[,assessment_year,1]-stf@m[,assessment_year,1]*stf@m.spwn[,assessment_year,1]) *
                                                     stf@mat[,assessment_year,1])),digits = 0),
                                RImY = rnd(drop(stf@stock.n[ac(0),assessment_year,1]),digits = 0),
                                RFcY = rnd(drop(stf@stock.n[ac(0),ac(an(assessment_year)+1),1]),digits = 0),
                                catchImY = rnd(drop(unitSums(stf@catch[,assessment_year])),digits = 0))

# ImY table
ImY_table <- forecast[forecast$Basis == 'intermediate year',]
ImY_table <- ImY_table[,ImY_columns]

# write outputs
write.taf(forecast_final, dir="report")
write.taf(ImY_assumptions, dir="report")
write.taf(ImY_table, dir="report")

# plot forecast options

filtOpt <- c('fmsyAR_notransfer',
             'fmsyAR_TACrule_transfer',
             'nf',
             'fmsy')

forecast_finalFilt <- forecast_final[forecast_final$Basis %in% filtOpt,]

taf.png("catch_forecast_options")
print(ggplot()+
        geom_bar(data=catchStf,aes(x=year,y=data),stat="identity")+
        geom_point(data=forecast_finalFilt,aes(x=an(assessment_year)+1,y=CatchA,col=Basis),size=2)+
        ylab('A fleet catches'))
dev.off()

taf.png("ssb_forecast_options")
print(ggplot()+
        geom_line(data=ssbStf,aes(x=year,y=data),linewidth=2)+
        geom_point(data=forecast_finalFilt,aes(x=an(assessment_year)+1,y=SSB1,col=Basis),size=2)+
        geom_point(data=forecast_finalFilt,aes(x=an(assessment_year)+2,y=SSB2,col=Basis),size=2)+
        geom_hline(yintercept = 1131601,col='blue')+
        geom_hline(yintercept = 903707,col='yellow')+
        geom_hline(yintercept = 828874,col='red')+
        ylab('NSAS SSB'))
dev.off()
