## Prepare plots for report

## Before: results_sf.RData (model)
## After:  summary.png (report)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scales)

mkdir("report")

data.source     <- file.path(".", "bootstrap",'data')
data.save       <- file.path(".", "data")
model.save      <- file.path(".", "model")
output.save     <- file.path(".", "output")

source('utilities_plots.R')

SMSkeyRuns  <- 2023
runName <- 'NSAS_HAWG2024'

load(file.path(model.save,paste0(runName,'_sf.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam
load(file.path(model.save,paste0(runName,'_sf_retro.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam
load(file.path(model.save,paste0(runName,'_mf.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam
load(file.path(model.save,paste0(runName,'_mf_retro.RData')))  # NSH, NSH.tun, NSH.ctrl, NSH.sam

# name(NSH) <- runName
# name(NSHs3$residual) <- runName
# name(NSHs3$sum) <- runName

################################################################################
### ============================================================================
### input data
### ============================================================================
################################################################################

prefix <- 'input'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))


##################################################
## stock weight at age
#################################################
taf.png(paste0(prefix,"_wstock_timeseries"))
zoom(timeseries(window(NSH,1975,range(NSH)["maxyear"]),slot="stock.wt"))
dev.off()

##################################################
## catch at age
#################################################
taf.png(paste0(prefix,"_wcatch_timeseries"))
zoom(timeseries(window(NSH,1975,range(NSH)["maxyear"]),slot="catch.wt"))
dev.off()

##################################################
## maturity at age
#################################################
taf.png(paste0(prefix,"_maturity"))
zoom(timeseries(window(NSH,1990,range(NSH)["maxyear"]),slot="mat"))
dev.off()

##################################################
## M at age
#################################################
taf.png(paste0(prefix,"_natmort"))
zoom(timeseries(window(NSH,1947,range(NSH)["maxyear"]),slot="m"))
dev.off()

##################################################
## prop index at age since 2000
#################################################
taf.png(paste0(prefix,"_prop_HERAS"))
prop_HERAS <- as.data.frame(NSH.tun[["HERAS"]]@index)
prop_HERAS$data[prop_HERAS$data == -1] <- NA

ggplot(subset(prop_HERAS,year>=2000),aes(x=year,y=data))+
  geom_bar(aes(fill = as.factor(age)),stat="identity",position = "fill")+
  ylab('Proportion of Acoustic index at age')+
  labs(fill = 'age (wr)')
dev.off()

##################################################
## HERAS ssb
#################################################
taf.png(paste0(prefix,"_ssb_HERAS_assessment"))
HERAS_SSB <- quantSums(NSH.tun[["HERAS"]]@index*
                        NSH@stock.wt[ac(1:8),ac(NSH.tun[["HERAS"]]@range['minyear']:NSH.tun[["HERAS"]]@range['maxyear'])]*
                        NSH@mat[ac(1:8),ac(NSH.tun[["HERAS"]]@range['minyear']:NSH.tun[["HERAS"]]@range['maxyear'])])
HERAS_SSB <- as.data.frame(HERAS_SSB)
HERAS_SSB$type <- 'HERAS'
assessment_ssb <- ssb(NSH.sam)
assessment_ssb$type <- 'assessment'
# ssb_df <- rbind(HERAS_SSB,assessment_ssb)

ggplot()+
  theme_bw()+
  geom_line(data=subset(HERAS_SSB,year>=2000),aes(x=year,y=data*1e-6,col='HERAS'),linewidth=1.5)+
  geom_line(data=subset(assessment_ssb,year>=2000),aes(x=year,y=value*1e-6,col='assessment'),linewidth=1.5)+
  geom_ribbon(data=subset(assessment_ssb,year>=2000),aes(x=year,ymin=lbnd*1e-6,ymax=ubnd*1e-6),alpha=0.5)+
  scale_colour_manual(values=c('grey','black'))+
  scale_fill_manual(values=c('grey'))+
  ylab('NSAS SSB (million t)')
dev.off()

##################################################
## prop catches at age since 2000
#################################################
taf.png(paste0(prefix,"_prop_catches"))
prop_catches <- as.data.frame(NSH@catch.n)
prop_catches$data[prop_catches$data == -1] <- NA

ggplot(subset(prop_catches,year>=2000),aes(x=year,y=data))+
  geom_bar(aes(fill = as.factor(age)),stat="identity",position = "fill")+
  ylab('Proportion of Catch numbers at age')+
  labs(fill = 'age (wr)')
dev.off()

##################################################
## weight in the stock by cohort
#################################################
taf.png(paste0(prefix,"_wstock_cohort"))
west.by.cohort      <- as.data.frame(FLCohort(window(NSH@stock.wt,2000,range(NSH)["maxyear"])))
west.by.cohort      <- subset(west.by.cohort,!is.na(west.by.cohort$data))
west.by.cohort$year <- west.by.cohort$age + west.by.cohort$cohort
zoom(xyplot(data~year,data=west.by.cohort,
            groups=cohort,
            auto.key=list(space="right",points=FALSE,lines=TRUE,type="b"),
            type="b",
            xlab="Year",ylab="Weight in the stock (kg)",
            main=paste(NSH@name,"Weight in the stock by cohort"),
            par.settings=list(superpose.symbol=list(pch=as.character(unique(west.by.cohort$cohort)%%10),cex=1.25)),
            panel=function(...) {
              panel.grid(h=-1,v=-1)
              panel.xyplot(...)
            }))
dev.off()

##################################################
## overlay survey time series
#################################################
start_year  <- 2000
end_year    <- dims(NSH)$maxyear+1

for(iAge in range(NSH)["min"]:range(NSH)["max"]){
  taf.png(paste0(prefix,"_overlay_survey_series_age_", iAge))
  for(idxSurvey in 1:length(NSH.tun)){
    # make sure we don't take any of the LAI data
    if(!grepl('LAI',name(NSH.tun[[idxSurvey]]))){
      idxFilt <- which(iAge == NSH.tun[[idxSurvey]]@range[1]:NSH.tun[[idxSurvey]]@range[2])
      
      if(idxSurvey == 1){
        dat <- cbind(as.data.frame(NSH.tun[[idxSurvey]]@index[idxFilt,]),
                     rep(NSH.tun[[idxSurvey]]@name,
                         length(drop(NSH.tun[[idxSurvey]]@index[idxFilt,]))))
        colnames(dat)[dim(dat)[2]] <- 'survey'
        dat$data <- (dat$data-mean(dat$data[dat$data != -1]))/sd(dat$data[dat$data != -1])
      }else{
        tempVar <- cbind(as.data.frame(NSH.tun[[idxSurvey]]@index[idxFilt,]),
                         rep(NSH.tun[[idxSurvey]]@name,
                             length(drop(NSH.tun[[idxSurvey]]@index[idxFilt,]))))
        colnames(tempVar)[dim(tempVar)[2]] <- 'survey'
        tempVar$data <- (tempVar$data-mean(tempVar$data[tempVar$data != -1]))/sd(tempVar$data[tempVar$data != -1])
        dat     <- rbind(dat,tempVar)
      }
    }
  }
  
  p <- ggplot(data = dat,mapping=aes(x=year,y=data,color=survey))+
    theme_bw() +
    geom_line()+
    xlim(start_year, end_year)+
    ylab('standardized index')+xlab('year')+
    ggtitle(paste('age ',iAge))
  
  print(p)
  dev.off()
}

##################################################
## IBTS0/IBTSQ1 relationship
#################################################
taf.png(paste0(prefix,"_IBTS0_vs_IBTSQ1"))

IBTSQ1      <- read.table(file.path('../../bootstrap/data','DATRAS_IBTSQ1.csv'), sep=",", header = TRUE) # read raw indices instead of standardized ones
IBTSQ1      <- IBTSQ1[IBTSQ1$IndexArea == "NS_Her",]
#IBTSQ1$Year <- IBTSQ1$Year-1

array_IBTSQ1 <- data.frame(cbind(IBTSQ1$Year, IBTSQ1$Age_1))
names(array_IBTSQ1) <- c("year","data")
array_IBTSQ1$yearclass <- array_IBTSQ1$year-2

array_IBTS0   <- subset(as.data.frame(NSH.tun), cname == "IBTS0" & slot == "index")
array_IBTS0$yearclass <- array_IBTS0$year-1

minYear <- max(c(min(array_IBTSQ1$yearclass), 
                 min(array_IBTS0$yearclass)))

array_IBTSQ1 <- array_IBTSQ1[array_IBTSQ1$yearclass >= minYear,]
array_IBTS0 <- array_IBTS0[array_IBTS0$yearclass >= minYear,]

#array_IBTSQ1 <- array_IBTSQ1[2:dim(array_IBTSQ1)[1],] # this because IBTS0 is the limiting factor

# construct the array to be plotted
tabMarkers <- data.frame(cbind(array_IBTSQ1$year,
                               array_IBTSQ1$data,
                               array_IBTS0[1:dim(array_IBTS0)[1]-1,]$data)) # exclude last year from the vector

names(tabMarkers) <- c("year", "IBTSQ1", "IBTS0")

ggplot(data=tabMarkers[1:dim(tabMarkers)[1]-1,], aes(x= IBTS0, y = IBTSQ1, label = year)) + 
  ylab("IBTSQ1 index") + 
  xlab("IBTS0 index") + 
  geom_text(colour = "red")+
  geom_vline(xintercept=array_IBTS0[dim(array_IBTS0)[1],]$data,
             colour="blue")+
  geom_text(aes(x=array_IBTS0[dim(array_IBTS0)[1],]$data, 
                label=array_IBTS0[dim(array_IBTS0)[1],]$year, 
                y=max(tabMarkers$IBTSQ1)- min(tabMarkers$IBTSQ1)/2+min(tabMarkers$IBTSQ1)), 
            colour="blue", 
            angle=90, 
            vjust = -0.3)

dev.off()

##################################################
## IBTS-Q3 combined
#################################################

load(file.path('../../bootstrap/initial/data/IBTSQ3/','result_Q3.RData'))
IBTSQ3.df <- res
IBTSQ3.df <- IBTSQ3.df %>% pivot_longer(!year,names_to = 'age',values_to = 'data')
IBTSQ3.df$age <- an(gsub("age","",IBTSQ3.df$age))
IBTSQ3.df <- subset(IBTSQ3.df,age >= 2)

IBTSQ3.ssb <- NSH@mat[ac(2:6),
              ac(min(IBTSQ3.df$year):max(an(IBTSQ3.df$year)))]*
  IBTSQ3.df %>% pivot_wider(names_from = 'year',values_from = 'data') %>% select(-c('age'))
IBTSQ3.ssb$age <- c(2:6)
IBTSQ3.ssb <- IBTSQ3.ssb %>% pivot_longer(!age,names_to ='year',values_to = 'data') %>% 
              group_by(year) %>% summarize(data=sum(data))
IBTSQ3.ssb$source <- 'IBTS-Q3'

load(file.path('../../bootstrap/initial/data/IBTSQ1/','result_Q1.RData'))
IBTSQ1.df <- res
IBTSQ1.df <- IBTSQ1.df %>% pivot_longer(!year,names_to = 'age',values_to = 'data')
IBTSQ1.df$age <- an(gsub("age","",IBTSQ1.df$age))
IBTSQ1.df <- subset(IBTSQ1.df,age >= 2)

IBTSQ1.ssb <- NSH@mat[ac(2:6),
                      ac(min(IBTSQ1.df$year):(max(IBTSQ1.df$year)-1))]*
  IBTSQ1.df %>% pivot_wider(names_from = 'year',values_from = 'data') %>% select(-c('age'))
IBTSQ1.ssb$age <- c(2:6)
IBTSQ1.ssb <- IBTSQ1.ssb %>% pivot_longer(!age,names_to ='year',values_to = 'data') %>% 
  group_by(year) %>% summarize(data=sum(data))
IBTSQ1.ssb$source <- 'IBTS-Q1'

HERAS.ssb <- as.data.frame(quantSums(NSH.tun$HERAS@index[ac(2:6),]*
                                        NSH@mat[ac(2:6),
                                                ac(min(an(dimnames(NSH.tun$HERAS)$year)):max(an(dimnames(NSH.tun$HERAS)$year)))]))
HERAS.ssb$source <- 'HERAS'
HERAS.ssb <- HERAS.ssb %>% select(colnames(IBTSQ1.ssb))

assessment.ssb <- as.data.frame(ssb(NSH[ac(2:6),]))
assessment.ssb$source <- 'assessment'
assessment.ssb <- assessment.ssb %>% select(colnames(IBTSQ1.ssb))


survey.all <- rbind(IBTSQ1.ssb,IBTSQ3.ssb,HERAS.ssb,assessment.ssb)
survey.all <- survey.all %>% group_by(source) %>% mutate(data = (data - mean(data))/sd(data))

taf.png(paste0(prefix,"_survey trends"))
ggplot(subset(survey.all,year >= 2000),aes(x=an(year),y=data,col=source))+
  geom_line()+
  ylab('z-normalized')+
  ggtitle('Normalized SSB trends (age 2-6)')
dev.off()

##################################################
## Internal consistency HERAS
#################################################
taf.png(paste0(prefix,"_HERAS_internal_consistency"))
zoom(plot(NSH.tun[["HERAS"]],type="internal"))
dev.off()

##################################################
## Internal consistency IBTSQ3
#################################################
taf.png(paste0(prefix,"_IBTSQ3_internal_consistency"))
zoom(plot(NSH.tun[["IBTS-Q3"]],type="internal"))
dev.off()

setwd('../..')

################################################################################
### ============================================================================
### Single-fleet
### ============================================================================
################################################################################

prefix <- 'sf'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

##################################################
## assessment fit
##################################################
pdf(file.path(paste0(prefix,"_fit.pdf")))
residual.diagnostics(NSH.sam)
dev.off()

##################################################
## stock trajectory
##################################################
taf.png(paste0(prefix,"_stock_trajectory"))
zoom(plot(NSH.sam))
dev.off()

##################################################
## Catchability at age
##################################################
taf.png(paste0(prefix,"_survey_catchability"))
catch <- catchabilities(NSH.sam)
zoom(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
            scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
            type=c("l","p"),lwd=c(2,1,1),col=c("black","grey","grey"),
            subset=fleet %in% c("HERAS","IBTS-Q1","IBTS-Q3"),
            main="Survey catchability parameters",ylab="Catchability",xlab="Age"))
dev.off()

##################################################
## observation variance by data source
##################################################
taf.png(paste0(prefix,"_observation_var_by_source"))
obv <- obs.var(NSH.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
dev.off()

##################################################
## variance vs uncertainty for each data source
#################################################
taf.png(paste0(prefix,"_variance_vs_uncertainty"))
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
dev.off()

##################################################
## Selectivity pattern per 5 years
#################################################
taf.png(paste0(prefix,"_selectivity"))
sel.pat <- merge(f(NSH.sam),fbar(NSH.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
zoom(xyplot(sel ~ age|sprintf("%i's",floor((year)/5)*5),sel.pat,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar"))
dev.off()

##################################################
## correlation matrix model params
#################################################
taf.png(paste0(prefix,"_cor_params"))
cor.plot(NSH.sam)
dev.off()

##################################################
## catch residuals per year per age
#################################################
taf.png(paste0(prefix,"_catage_residuals"))
dat <- subset(residuals(NSH.sam),fleet=="catch unique")
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=1*abs(lst$cex))
            }))
dev.off()

##################################################
## IBTS-Q1 residuals per year per age
#################################################
taf.png(paste0(prefix,"_IBTSQ1_residuals"))
dat <- subset(residuals(NSH.sam),fleet=="IBTS-Q1")
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year IBTSQ1",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=1*abs(lst$cex))
            }))
dev.off()

##################################################
## HERAS residuals per year per age
#################################################
taf.png(paste0(prefix,"_HERAS_residuals"))
dat <- subset(residuals(NSH.sam),fleet=="HERAS")
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year HERAS",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=1*abs(lst$cex))
            }))
dev.off()

##################################################
## IBTS-Q3 residuals per year per age
#################################################
taf.png(paste0(prefix,"_IBTSQ3_residuals"))
dat <- subset(residuals(NSH.sam),fleet=="IBTS-Q3")
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year IBTS-Q3",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=1*abs(lst$cex))
            }))
dev.off()

##################################################
## proc error N
#################################################
taf.png(paste0(prefix,"_procErr_N"))
procerr.plot(NSH+NSH.sam,weight="stock.wt",type="n",rel=T)
dev.off()

##################################################
## proc error M
#################################################
taf.png(paste0(prefix,"_procErr_M"))
procerr.plot(NSH+NSH.sam,weight="stock.wt",type="mort",rel=T)
dev.off()


##################################################
## F at age
#################################################
taf.png(paste0(prefix,"_fatage"))
zoom(timeseries(window(NSH,2000,range(NSH)["maxyear"]),slot="harvest"))
dev.off()

##################################################
## F/fbar at age
#################################################
NSH.temp <- NSH
NSH.temp@harvest <- sweep(NSH@harvest,2:6,fbar(NSH),'/')

taf.png(paste0(prefix,"_ffbaratage"))
zoom(timeseries(window(NSH.temp,2000,range(NSH)["maxyear"]),slot="harvest"))
dev.off()

##################################################
## Spawning component proportions
#################################################
taf.png(paste0(prefix,"_component_proportions"))
components.df <- as.data.frame(components(NSH.sam))
components.df <- components.df[,1:dim(components.df)[2]]
colnames(components.df) <- ac(1972:range(NSH)['maxyear'])
ssb.df <- ssb(NSH.sam)
ssb.df <- ssb.df[ssb.df$year %in% 1972:range(NSH)['maxyear'],]
components_ssb.df <- sweep(components.df, 2, ssb.df$value,'*')
components.df$components <- rownames(components.df)
components.df <- components.df %>% pivot_longer(!components,names_to = "year",values_to='Fraction')
components.df$year <- as.numeric(components.df$year)
ggplot(components.df,aes(x=year,y=Fraction,fill=components))+
  geom_area()
dev.off()

##################################################
## Spawning component SSB
#################################################

taf.png(paste0(prefix,"_component_SSB"))
components_ssb.df['all',] <- apply(components_ssb.df,2,'sum')
components_ssb.df$components <- rownames(components_ssb.df)
components_ssb.df <- components_ssb.df %>% pivot_longer(!components,names_to = "year",values_to='ssb')
components_ssb.df$year <- as.numeric(components_ssb.df$year)
ggplot(components_ssb.df,aes(x=year,y=ssb,col=components))+
  geom_line()
dev.off()

##################################################
## Recruitment vs SSB
#################################################

df.rec <- rec(NSH.sam)
df.rec <- df.rec %>% select(c('year','value'))
colnames(df.rec) <- c('year','rec')
df.ssb <- ssb(NSH.sam)
df.ssb <- df.ssb %>% select(c('year','value'))
colnames(df.ssb) <- c('year','ssb')

df.all <- left_join(df.rec,df.ssb,by=('year'))
df.all$yearclass <- df.all$year-1
df.all$ssbR <- df.all$ssb*100/df.all$rec
df.all$log.ssbR <- log10(df.all$ssbR)
df.all$RperS <- df.all$rec/df.all$ssb

df.IBTS0 <- as.data.frame(NSH.tun$IBTS0@index)
df.IBTS0 <- subset(df.IBTS0,year<=(max(df.IBTS0$year)-1))
df.IBTS0 <- df.IBTS0 %>% select(c('year','data'))
colnames(df.IBTS0) <- c('year','IBTS0')
df.IBTS0 <- left_join(df.IBTS0,subset(df.ssb,year %in% df.IBTS0$year),by=('year'))
df.IBTS0$larSurTot <- df.IBTS0$IBTS0/(df.IBTS0$ssb/1000)

prop.comp <- as.data.frame(components(NSH.sam))
prop.comp$components <- rownames(prop.comp)
prop.comp <- prop.comp %>% 
  pivot_longer(!components,names_to = 'year',values_to = 'prop')%>% 
  pivot_wider(names_from = components,values_from = prop)

prop.comp$year <- as.numeric(prop.comp$year)

df.IBTS0 <- left_join(df.IBTS0,subset(prop.comp,year %in% df.IBTS0$year),by=('year'))
df.IBTS0$larSurNorth <- df.IBTS0$IBTS0/(df.IBTS0$ssb/1000)*(df.IBTS0$`LAI-ORSH`+df.IBTS0$`LAI-BUN`+df.IBTS0$`LAI-CNS`)
df.IBTS0$larSurSouth <- df.IBTS0$IBTS0/(df.IBTS0$ssb/1000)*(df.IBTS0$`LAI-SNS`)
df.IBTS0$yearclass <- df.IBTS0$year-1
df.IBTS0 <- df.IBTS0 %>% select(c(yearclass,larSurSouth,larSurNorth,larSurTot)) %>%
  pivot_longer(!yearclass,names_to = 'Component',values_to = 'larSur')
df.IBTS0$Component[df.IBTS0$Component == "larSurSouth"] <- 'Southern'
df.IBTS0$Component[df.IBTS0$Component == "larSurNorth"] <- 'Northern'
df.IBTS0$Component[df.IBTS0$Component == "larSurTot"] <- 'Total'

taf.png(paste0(prefix,"_SR_periods"))
ggplot(data=df.all,aes(x=ssb,y=rec,col=year>2001),aes(x=ssb,y=rec,col=year>2001))+
  geom_point(size=3)+
  geom_point(data=subset(df.all,year == max(df.all$year)),aes(x=ssb,y=rec),col='grey20', size=3)+
  scale_color_manual(values = c("red","blue"))+
  ylab('Recruitment 0-wr (thousands)')+
  xlab('Spawning Stock Biomass (t)')+
  theme_bw()
dev.off()

taf.png(paste0(prefix,"_recruit_spawners"))
ggplot()+
  geom_smooth(data=df.all,aes(x=yearclass,y=RperS),span = 0.1,col='black', linewidth=2)+
  geom_point(data=df.all,aes(x=yearclass,y=RperS,col=year>2001), size=2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()+
  annotation_logticks()+
  scale_color_manual(values = c("red","blue"))+
  xlab('Year class')+
  ylab('Recruits per spawner')
dev.off()

taf.png(paste0(prefix,"_larval_survival"))
df.IBTS0 %>% filter(Component != "Southern") %>% 
  ggplot(aes(x=yearclass,y=larSur,col=Component))+
  geom_point(size=3)+
  geom_line(linewidth=2)+
  #geom_smooth(span = 0.2,se=F)+
  theme_bw()+
  scale_color_manual(values = c("red","blue"))+
  xlab('Year class')+
  ylab('Larval survival')
dev.off() 

##################################################
## retro stock trajectory
#################################################
taf.png(paste0(prefix,"_retrospective_stock"))
zoom(plot(NSH.retro))
dev.off()

##################################################
## model parameters retrospective
#################################################
taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
zoom(retroParams(NSH.retro))
dev.off()

##################################################
## fishing selectivity retrospective
#################################################
taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
zoom(retroSelectivity(NSH.retro,(range(NSH)["maxyear"]-10):range(NSH)["maxyear"]))
dev.off()

##################################################
## FMSYAR
#################################################
taf.png(paste0(prefix,"_FMSYAR"))
plot(x=c(0,1.23,2.6),
     y=c(0,0.3,0.3),
     type="l",
     ylim=c(0,0.4),
     lwd=2,
     xlab="SSB in million tonnes",
     ylab="Fbar",
     cex.lab=1.3,
     main="FMSYAR North Sea Herring")
abline(v=0.87,col="red",lwd=2,lty=2)
abline(v=0.96,col="blue",lwd=2,lty=2)
abline(v=1.23,col="darkgreen",lwd=2,lty=2)
text(0.87,0,labels=expression(B[lim]),col="red",cex=1.3,pos=2)
text(0.96,0,labels=expression(B[pa]),col="blue",cex=1.3,pos=2)
text(1.23,0,labels=expression(B[trigger]),col="darkgreen",cex=1.3,pos=4)

points(y=fbar(NSH[ac(2:6),ac((range(NSH)[5]-10):range(NSH)[5])]), x=(ssb(NSH[,ac((range(NSH)[5]-10):range(NSH)[5])])/1e6),pch=19)
lines(y=fbar(NSH[ac(2:6),ac((range(NSH)[5]-10):range(NSH)[5])]),  x=(ssb(NSH[,ac((range(NSH)[5]-10):range(NSH)[5])])/1e6))
text(y=fbar(NSH[ac(2:6),ac((range(NSH)[5]-10):range(NSH)[5])]),   x=(ssb(NSH[,ac((range(NSH)[5]-10):range(NSH)[5])])/1e6),
     labels=ac((range(NSH)[5]-10):range(NSH)[5]),pos=3,cex=0.7)
dev.off()



setwd('../..')

################################################################################
### ============================================================================
### Multi-fleet
### ============================================================================
################################################################################

prefix <- 'mf'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))


##################################################
## stock trajectory
##################################################
taf.png(paste0(prefix,"_stock_trajectory"))
plot(NSH3f.sam)
dev.off()

##################################################
## catchability at age mf
##################################################
taf.png(paste0(prefix,"_survey_catchability"))
catch <- catchabilities(NSH3f.sam)
xyplot(value+ubnd+lbnd ~ age | fleet,catch,
       scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
       type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
       subset=fleet %in% c("HERAS","IBTS-Q3"),
       main="Survey catchability parameters",ylab="Catchability",xlab="Age")
dev.off()

##################################################
## variance per data source
##################################################
taf.png(paste0(prefix,"_observation_var_by_source"))
obv <- obs.var(NSH3f.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
dev.off()

##################################################
## obs.var vs uncertainty mf
##################################################
taf.png(paste0(prefix,"_variance_vs_uncertainty"))
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
dev.off()

##################################################
## selectivity A fleet
##################################################
sel.pat <- as.data.frame(sweep(NSH3f.sam@harvest[,,,,'A'],2:6,quantMeans(NSH3f.sam@harvest[ac(2:6),,,,'A']),"/"))
sel.pat <- rbind(sel.pat,
                 as.data.frame(sweep(NSH3f.sam@harvest[,,,,'BD'],2:6,quantMeans(NSH3f.sam@harvest[ac(0:1),,,,'BD']),"/")))
sel.pat <- rbind(sel.pat,
                 as.data.frame(sweep(NSH3f.sam@harvest[,,,,'C'],2:6,quantMeans(NSH3f.sam@harvest[ac(1:3),,,,'C']),"/")))

taf.png(paste0(prefix,"_selectivity_A"))

sel.ind <- subset(sel.pat,area=='A')
sel.ind$sel <- sel.ind$data
xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.ind,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="A fleet - Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar")
dev.off()

##################################################
## selectivity C fleet
##################################################
taf.png(paste0(prefix,"_selectivity_C"))

sel.ind <- subset(sel.pat,area=='C')
sel.ind$sel <- sel.ind$data
xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.ind,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="C fleet - Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar")
dev.off()

##################################################
## selectivity BD fleet
##################################################

taf.png(paste0(prefix,"_selectivity_BD"))
sel.ind <- subset(sel.pat,area=='BD')
sel.ind$sel <- sel.ind$data
xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.ind,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="BD fleet - Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar")
dev.off()

##################################################
## recent selectivity A fleet
##################################################
sel.pat <- merge(f(NSH3f.sam),fbar(NSH3f.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))

taf.png(paste0(prefix,"_selectivity_recent_A"))
# Fleet A
selA <- subset(sel.pat, fleet == "catch A")
selA <- subset(selA, year > 2012)
for(idxYear in 1:length(selA$year)){selA$year[idxYear] <- toString(selA$year[idxYear])}

print(ggplot(selA, aes(x = age, y = sel, colour = year)) + geom_line() + xlab("age") + ylab("F/Fbar") + ggtitle("fleet A"))

dev.off()

##################################################
## recent selectivity C fleet
##################################################
taf.png(paste0(prefix,"_selectivity_recent_C"))
selC <- subset(sel.pat, fleet == "catch C")
selC <- subset(selC, year > 2012)
for(idxYear in 1:length(selC$year)){selC$year[idxYear] <- toString(selC$year[idxYear])}

print(ggplot(selC, aes(x = age, y = sel, colour = year)) + geom_line() + xlab("age") + ylab("F/Fbar") + ggtitle("fleet C"))
dev.off()

##################################################
## recent selectivity BD fleet
##################################################
taf.png(paste0(prefix,"_selectivity_recent_BD"))
selBD <- subset(sel.pat, fleet == "catch BD")
selBD <- subset(selBD, year > 2012 & age <= 2)
for(idxYear in 1:length(selBD$year)){selBD$year[idxYear] <- toString(selBD$year[idxYear])}

print(ggplot(selBD, aes(x = age, y = sel, colour = year)) + geom_line() + xlab("age") + ylab("F/Fbar") + ggtitle("fleet BD"))
dev.off()

##################################################
## parameters correlation
##################################################
taf.png(paste0(prefix,"_cor_params"))
cor.plot(NSH.sam)
dev.off()

##################################################
## catage_residuals A residuals
##################################################
taf.png(paste0(prefix,"_catage_residuals_A"))
dat <- subset(residuals(NSH3f.sam),fleet=="catch A")
dat <- subset(dat, year >= 1997)
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch fleet A",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=3*abs(lst$cex))
            }))
dev.off()

##################################################
## catage_residuals C residuals
##################################################
taf.png(paste0(prefix,"_catage_residuals_C"))
dat <- subset(residuals(NSH3f.sam),fleet=="catch C")
dat <- subset(dat, year >= 1997)
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch fleet C",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=3*abs(lst$cex))
            }))
dev.off()

##################################################
## catage_residuals BD residuals
##################################################
taf.png(paste0(prefix,"_catage_residuals_BD"))
dat <- subset(residuals(NSH3f.sam),fleet=="catch BD")
dat <- subset(dat, year >= 1997 & age <= 2)
zoom(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch fleet BD",
            panel=function(...){
              lst <- list(...)
              panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=3*abs(lst$cex))
            }))
dev.off()

##################################################
## sf/mf comparison
##################################################
NSH.sams  <- new("FLSAMs")
NSH.sams[['sf']] <- NSH.sam
NSH.sams[['mf']] <- NSH3f.sam

out_assessment <- rbind(cbind(ssb(NSH.sams),type='ssb'),
                        cbind(rec(NSH.sams),type='rec'),
                        cbind(fbar(NSH.sams),type='fbar'))

taf.png(paste0(prefix,"_comp_models recent"))
print(ggplot(subset(out_assessment,year>=2000),aes(x=year,y=value,fill=name,col=name))+
        geom_line()+
        geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
        facet_wrap(~type,scales='free',ncol=1))
dev.off()

taf.png(paste0(prefix,"_comp_models full"))
print(ggplot(out_assessment,aes(x=year,y=value,fill=name,col=name))+
        geom_line()+
        geom_ribbon(aes(ymin=lbnd,ymax=ubnd),alpha=0.5, colour = NA)+
        facet_wrap(~type,scales='free',ncol=1))
dev.off()

##################################################
## retro stock trajectory
#################################################
taf.png(paste0(prefix,"_retrospective_stock"))
zoom(plot(NSH3f.retro))
dev.off()

setwd('../..')
