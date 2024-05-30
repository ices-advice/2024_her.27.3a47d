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

mkdir("report")

data.source <- file.path('.','bootstrap','data')
model.dir  <- file.path('.','model')

sourceDir("bootstrap/initial/software/utilities")


## -----------------------------------------------------------------------------
## plotting individual stfs
## -----------------------------------------------------------------------------
assessment_year     <- '2024'

# load forecast results (stf & stf2.df)
load(file.path(model.dir,paste0('forecast_',assessment_year,'.RData')))
stf2.df <- stf2.df %>% mutate(source="forecast")

# load assessment (NSH, NSH.ctrl, NSH.sam, NSH.tun)
load(file.path(data.source,paste0('NSAS_HAWG',assessment_year,'_sf.RData')))

# load TAC
TACTab            <- read.taf("bootstrap/data/stf_NSAS_TAC.csv")
adviceTab         <- read.taf("bootstrap/data/stf_NSAS_advice.csv")

adviceTab <- pivot_longer(adviceTab,!year,names_to='fleets',values_to='data')
adviceTab$type <- 'advise'
TACTab <- pivot_longer(TACTab,!year,names_to='fleets',values_to='data')
TACTab$type <- 'TAC'
adviseTAC <- rbind(adviceTab,TACTab)

# try(setwd(path),silent=FALSE)

# get years for previous and current forecasts
y1 <- an(range(stf.res[1,2]$stf[[1]])[["minyear"]])+1
y2 <- y1+1
y3 <- y2+1

stf_y3 <- stf.res[stf.res$scenario=="fmsyAR_notransfer",2]$stf[[1]]

stf_plot_names    <- paste0(assessment_year,'_stf')

# load previous year forecasts
stf_y1 <- extract_from_rdata(paste0("bootstrap/data/forecast_",y1,".RData"), "stf.res")
stf_y1 <- stf_y1[stf_y1$scenario=="fmsyAR_no_transfer",2]$stf[[1]]
stf_y2 <- extract_from_rdata(paste0("bootstrap/data/forecast_",y2,".RData"), "stf.res")
stf_y2 <- stf_y2[stf_y2$scenario=="fmsyAR_no_transfer",2]$stf[[1]]

# set wd to report
setwd("report")

#Read in data (this is not really needed)
eval(expr=parse(text=paste0('resy1 <- stf_y1')))
eval(expr=parse(text=paste0('resy2 <- stf_y2')))
eval(expr=parse(text=paste0('resy3 <- stf_y3')))

resy1.df <- as.data.frame(resy1) %>% mutate(stf=ac(y1))
resy2.df <- as.data.frame(resy2) %>% mutate(stf=ac(y2))
resy3.df <- as.data.frame(resy3) %>% mutate(stf=ac(y3))

# THIS IS NOT THE RIGHT FORMULA !! Pick up results from df instead
# ssby1    <- as.data.frame(resy1@stock.n * resy1@stock.wt * resy1@mat) %>% mutate(stf=ac(y1))
# ssby2    <- as.data.frame(resy2@stock.n * resy2@stock.wt * resy2@mat) %>% mutate(stf=ac(y2))
# ssby3    <- as.data.frame(resy3@stock.n * resy3@stock.wt * resy3@mat) %>% mutate(stf=ac(y3))

df.catches.y2 <- as.data.frame(stf_y2@catch[,ac(2022)])
df.catches.y2$assessment <- y2
df.catches.y3 <- as.data.frame(stf_y3@catch[,ac(2022)])
df.catches.y3$assessment <- y3

df.catches <- rbind(df.catches.y2,df.catches.y3)
df.catches <- df.catches %>% pivot_wider(values_from = data,names_from = assessment)

# windows()
# ggplot(df.catches,aes(x=year,y=data,col=as.factor(assessment)))+
#   geom_line()+
#   facet_wrap(~unit)

ssby1    <- 
  as.data.frame(resy1@stock.n[,,1] * resy1@stock.wt[,,1] *
                exp(-unitSums(resy1@harvest[,])*resy1@harvest.spwn[,,1]-resy1@m[,,1] *
                resy1@m.spwn[,,1]) * resy1@mat[,,1]) %>% 
  mutate(stf=ac(y1))
ssby2    <- 
  as.data.frame(resy2@stock.n[,,1] * resy2@stock.wt[,,1] *
                  exp(-unitSums(resy2@harvest[,])*resy2@harvest.spwn[,,1]-resy2@m[,,1] *
                        resy2@m.spwn[,,1]) * resy2@mat[,,1]) %>% 
  mutate(stf=ac(y2))
ssby3    <- 
  as.data.frame(resy3@stock.n[,,1] * resy3@stock.wt[,,1] *
                  exp(-unitSums(resy3@harvest[,])*resy3@harvest.spwn[,,1]-resy3@m[,,1] *
                        resy3@m.spwn[,,1]) * resy3@mat[,,1]) %>% 
  mutate(stf=ac(y3))

stocky1  <- as.data.frame(resy1@stock.n * resy1@stock.wt) %>% mutate(stf=ac(y1))
stocky2  <- as.data.frame(resy2@stock.n * resy2@stock.wt) %>% mutate(stf=ac(y2))
stocky3  <- as.data.frame(resy3@stock.n * resy3@stock.wt) %>% mutate(stf=ac(y3))

#catchy1  <- as.data.frame(resy1@catch.n * resy1@catch.wt) %>% mutate(stf=ac(y1)) %>% filter(year >= as.numeric(stf), !is.na(data))
#catchy2  <- as.data.frame(resy2@catch.n * resy2@catch.wt) %>% mutate(stf=ac(y2)) %>% filter(year >= as.numeric(stf), !is.na(data))
#catchy3  <- as.data.frame(resy3@catch.n * resy3@catch.wt) %>% mutate(stf=ac(y3)) %>% filter(year >= as.numeric(stf), !is.na(data))
catchy1    <- 
  as.data.frame(resy1@catch.n[,,1] * resy1@catch.wt[,,1]) %>% 
  mutate(stf=ac(y1))
catchy2    <- 
  as.data.frame(resy2@catch.n[,,1] * resy2@catch.wt[,,1]) %>% 
  mutate(stf=ac(y2))
catchy3    <- 
  as.data.frame(resy3@catch.n[,,1] * resy3@catch.wt[,,1]) %>% 
  mutate(stf=ac(y3))

harvesty1  <- as.data.frame(resy1@harvest) %>% mutate(stf=ac(y1))
harvesty2  <- as.data.frame(resy2@harvest) %>% mutate(stf=ac(y2))
harvesty3  <- as.data.frame(resy3@harvest) %>% mutate(stf=ac(y3))

catchAssessment <- as.data.frame(NSH@catch)

NSH.df <- 
  as.data.frame(NSH) %>% 
  bind_rows(as.data.frame(ssb(NSH)) %>% mutate(slot="ssb")) %>% 
  bind_rows(as.data.frame(rec(NSH)) %>% mutate(slot="rec", age=as.character(age))) %>% 
  bind_rows(as.data.frame(fbar(NSH)) %>% mutate(slot="fbar")) %>% 
  mutate(source="assessment")

##################################################
# plot SSB at age - lines
##################################################
taf.png("ssb_at_age")

print(bind_rows(ssby1,ssby2,ssby3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_line(aes(colour=stf), size=1) +
        labs(title="SSB at age", y="SSB (tonnes)") +
        facet_wrap(~age))

dev.off()

taf.png("ssb_at_age_bars")
print(bind_rows(ssby1,ssby2,ssby3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity") +
        labs(title="SSB at age", fill="age (WR)", y="SSB (tonnes)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

##################################################
# plot stockbiomass at age - lines
##################################################
taf.png("stock_B_at_age")
print(bind_rows(stocky1,stocky2,stocky3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_line(aes(colour=stf), size=1) +
        labs(title="TSB at age", y="TSB (tonnes)") +
        facet_wrap(~age))
dev.off()

# plot stockbiomass at age - bars
taf.png("stock_B_bars")
print(bind_rows(stocky1,stocky2,stocky3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity") +
        labs(title="TSB at age", fill="age (WR)", y="TSB (tonnes)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

# plot stockbiomass at age - bars
taf.png("stock_B_bars_bySTF")
print(bind_rows(stocky1,stocky2,stocky3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=stf, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity") +
        labs(title="TSB at age by STF", fill="age (WR)", y="TSB (tonnes)") +
        facet_wrap(~year))

dev.off()

# plot catch and TAC
taf.png("catch_TAC")
print(ggplot()+
              geom_bar(data=catchAssessment,aes(x=year,y=data),stat="identity")+
              geom_line(data=subset(adviseTAC,fleets=='A'),aes(x=year,y=data,col=type),size=1.5)+
              xlim(NSH@range['maxyear']-15,NSH@range['maxyear']+1))
dev.off()

# plot catch and TAC
taf.png("TAC_advise")
print(ggplot(subset(adviseTAC,year>=2015 & fleets != 'F'),aes(x=year,y=data,col=type))+
                     geom_line()+
              facet_wrap(~fleets,scales='free'))
dev.off()

# plot harvest at age - line
taf.png("harvest_trajectory")
print(bind_rows(resy1.df,resy2.df, resy3.df) %>%
        filter(unit=="A") %>% 
        filter(slot=="harvest") %>% 
        ggplot(aes(x=year, y=data)) +
        theme_bw() +
        geom_line(aes(colour=factor(age))) +
        labs(title="harvest", fill="age (WR)") +
        facet_wrap(~paste("STF", stf)))
dev.off()

# plot catch biomass at age - line
taf.png("catch_age")
print(bind_rows(catchy1,catchy2,catchy3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_line(aes(colour=stf), size=1) +
        labs(title="Catch at age (tonnes)", y= "Catch (tonnes)", fill="STF") +
        facet_wrap(~paste(age, "WR")))

dev.off()

# plot catch biomass at age - bar
taf.png("catch_bars")
print(bind_rows(catchy1, catchy2,catchy3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity", colour="white") +
        labs(title="Catch at age (tonnes)", y = "Catch (tonnes)", fill="age (wr)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

# plot catch biomass at age - percentage bar
taf.png("catch_proportions")
print(bind_rows(catchy1,catchy2,catchy3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity", position="fill", colour="white") +
        scale_y_continuous(labels = scales::percent) +
        labs(title="Catch age age (tonnes, as proportions)", y="proportion", fill="age (wr)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

# plot catch biomass at age by fleet
taf.png("catch_bars_byfleet")
print(bind_rows(catchy1, catchy2,catchy3) %>% 
        # filter(as.numeric(year) >= as.numeric(stf)) %>% 
        # filter(unit=="B") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_bar(aes(fill=as.character(age)), stat="identity", colour="white") +
        labs(title="Catch at age (tonnes)", y = "Catch (tonnes)", fill="age (wr)") +
        facet_grid(unit~paste("STF", stf), scales="free_y"))

dev.off()

# catch weight age age - line
taf.png("catch_weight_A_fleet")
print(bind_rows(mutate(as.data.frame(resy1@catch.wt), stf=ac(y1)),
                mutate(as.data.frame(resy2@catch.wt), stf=ac(y2)),
                mutate(as.data.frame(resy3@catch.wt), stf=ac(y3)) ) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_line(aes(group=age, colour=as.character(age)), size=1) +
        labs(title="catch weight at age (kg, A fleet)",y="kg", colour="age (wr)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

taf.png("stock_weight")
print(bind_rows(mutate(as.data.frame(resy1@stock.wt[,,1]), stf=ac(y1)),
                mutate(as.data.frame(resy2@stock.wt[,,1]), stf=ac(y2)),
                mutate(as.data.frame(resy3@stock.wt[,,1]), stf=ac(y3)) ) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        geom_line(aes(group=age, colour=as.character(age)), size=1) +
        labs(title="Stock weight at age (kg)",y="kg", colour="age (wr)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

# catch number at age - line
taf.png("catch_trajectory")
print(bind_rows(mutate(as.data.frame(resy1@catch.n), stf=ac(y1)), 
                mutate(as.data.frame(resy2@catch.n), stf=ac(y2)), 
                mutate(as.data.frame(resy3@catch.n), stf=ac(y3))) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=stf)) +
        theme_bw() +
        expand_limits(y=0) +
        geom_line(aes(colour=stf), size=1) +
        labs(title="Catch numbers at age (thousands)",y="thousands", colour="age (wr)") +
        facet_wrap(~paste(age, "wr"), scales="free_y"))

dev.off()

# plot harvest at age - lines
taf.png("stock_harvest_at_age")
print(bind_rows(harvesty1,harvesty2,harvesty3) %>% 
        filter(unit=="A") %>% 
        ggplot(aes(x=year, y=data, group=age)) +
        theme_bw() +
        geom_line(aes(colour=factor(age)), size=1) +
        labs(title="F at age",y="F", colour="age (wr)") +
        facet_wrap(~paste("STF", stf)))

dev.off()

# plot medium term forecasts (without and with assessment)
taf.png("stf_MSYAR_projection")
print(stf2.df %>% 
        filter(slot %in% c("stock","catch", "FA2-6", "rec"), 
               unit=="A") %>%
        
        ggplot(aes(x=year,y=data)) +
        theme_bw() +
        theme(legend.position = "none") +
        geom_line(aes(colour=slot), size=1.3) +
        scale_x_continuous(breaks=unique(stf2.df$year)) +
        expand_limits(y=0) +
        labs(title="Medium term forecast MSY AR without transfer", y="",x="") +
        facet_wrap(~slot, scales="free_y"))
dev.off()

taf.png("stf_MSYAR_projection2")
print(bind_rows(stf2.df, NSH.df) %>% 
      filter(slot %in% c("stock","catch", "FA2-6","fbar","ssb", "rec"), 
             unit %in% c("A", "unique"), 
             year >= as.numeric(assessment_year) - 10) %>%
      mutate(slot = ifelse(slot=="FA2-6", "fbar", slot),
             slot = ifelse(slot=="stock","ssb", slot)) %>% 
      drop_na() %>% 
      
      ggplot(aes(x=year,y=data)) +
      theme_bw() +
      theme(legend.position = "none") +
      geom_rect(aes(xmin=as.numeric(assessment_year), xmax=max(stf2.df$year), ymin=0, ymax=Inf), fill="lightgray", alpha=0.4) +
      geom_line(aes(colour=slot, size=source, linetype=source)) +
      scale_x_continuous(breaks=scales::pretty_breaks()) +
      scale_size_manual(values=c(1.0, 1.3)) +
      scale_linetype_manual(values=c("dashed", "solid")) +
      expand_limits(y=0) +
      labs(title="Assessment and medium term forecast MSY AR without transfer", y="",x="") +
      facet_wrap(~slot, scales="free_y")
)

dev.off()

df.ssb        <- as.data.frame(ssb_multiFleet(stf_y3))
df.ssb$type   <- 'ssb'
df.fbar       <- as.data.frame(unitSums(fbar(stf_y3)))
df.fbar$type  <- 'fbar'
df.all       <- rbind(df.ssb,df.fbar)
df.all$assess_year <- assessment_year
df.all$yearType <- c('oY','oY','DtY','ImY','FcY','CtY')

df.ssb        <- as.data.frame(ssb_multiFleet(stf_y2))
df.ssb$type   <- 'ssb'
df.fbar       <- as.data.frame(unitSums(fbar(stf_y2)))
df.fbar$type  <- 'fbar'
df.all_y2       <- rbind(df.ssb,df.fbar)
df.all_y2$assess_year <- ac(an(assessment_year)-1)
df.all_y2$yearType <- c('oY','oY','DtY','ImY','FcY','CtY')

df.all <- rbind(df.all,df.all_y2)

df.all <- df.all %>% select(-unit) %>% pivot_wider(names_from = type, values_from = data)

# ssb intermediate year assumptions
taf.png("SSB_comp")
print(ggplot(df.all,aes(x=year,y=ssb*1e-6,col=assess_year))+
        geom_line()+
        geom_point(aes(shape=yearType))+
        ylab('SSB (million tonnes)'))
dev.off()


# FMSYAR function

refptsAll <- list(Fmsy = 0.32,
                  Fsq  = NA,
                  Flim = 0.39,
                  Fpa  = 0.33,
                  Blim = 828874,
                  Bpa  = 903707,
                  MSYBtrigger = 1131601,
                  Ftarget  = 0, # set this for individual cases
                  F01   = 0.05,
                  Btrigger = 0)

df.refpts <- data.frame(matrix(ncol = 2, nrow = 0))

df.refpts           <- rbind(df.refpts,cbind(c(0,refptsAll$MSYBtrigger*1e-6,2.6),c(0,refptsAll$Fmsy,refptsAll$Fmsy)))
colnames(df.refpts) <- c('ssb','FMSY')
df.refpts$ssb       <- as.numeric(df.refpts$ssb)
df.refpts$FMSY      <- as.numeric(df.refpts$FMSY)

taf.png("FMSYAR")
print(ggplot()+
        geom_line(data=df.refpts,aes(x=ssb,y=FMSY))+
        geom_text(data=subset(df.all,yearType!='CtY' & assess_year == assessment_year),aes(x=ssb*1e-6,y=fbar,label=year),hjust=0, vjust=0, col="grey")+
        geom_vline(xintercept = refptsAll$MSYBtrigger/1000000, col="green")+
        geom_vline(xintercept = refptsAll$Bpa/1000000, col="blue")+
        geom_vline(xintercept = refptsAll$Blim/1000000, col="red")+
        xlab('SSB (millions tonnes)')+
        theme_bw())
dev.off()

setwd("..")

