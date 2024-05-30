## Prepare plots for report

## Before: results_sf.RData (model)
## After:  summary.png (report)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(doBy)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(gridExtra)

mkdir("report")

data.source     <- file.path(".", "bootstrap",'data')
data.save       <- file.path(".", "data")
model.save      <- file.path(".", "model")
output.save     <- file.path(".", "output")

source('utilities_model.R')

SMSkeyRunsMat <- c(2010,2013,2016,2019,2023)


### ============================================================================
### compare M vectors
### ============================================================================
prefix <- 'new_natmort'

mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

M.raw <- read.csv(file.path('../../bootstrap/data','SMS_NSAS_M_raw.csv'))

taf.png(file.path(paste0(prefix,"_znorm raw1")))
M.raw %>%
      group_by(Source, age) %>% 
      mutate(M_std = (M - mean(M))/sd(M) ) %>% 
      ggplot(aes(x=year, y=M_std, group=as.character(age))) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(aes(colour=as.character(age))) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(colour = "Age (WR)") +
      guides(colour = guide_legend(nrow = 1)) +
      facet_wrap(~Source, ncol=1)
dev.off()

taf.png(file.path(paste0(prefix,"_znorm raw2")))
M.raw %>%
      group_by(Source, age) %>% 
      mutate(M_std = (M - mean(M))/sd(M) ) %>% 
      
      ggplot(aes(x=year, y=M_std, group=as.character(Source))) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(aes(colour=as.character(Source))) +
      # scale_colour_manual(values=myPlotColors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(colour = "Keyrun") +
      guides(colour = guide_legend(nrow = 1)) +
      facet_wrap(~age, nrow=3)
dev.off()

taf.png(file.path(paste0(prefix,"_M raw")))
# M by Key run (facet) and WR (colour)
M.raw %>%
      ggplot(aes(x=year, y=M, group=as.character(age))) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(aes(colour=as.character(age))) +
      labs(colour = "Age (WR)") +
      # scale_colour_manual(values=myPlotColors) +
      guides(colour = guide_legend(nrow = 1)) +
      facet_wrap(~Source, nrow=1)
dev.off()

taf.png(file.path(paste0(prefix,"_M raw2")))
# M by WR (facet) and Key run (colour)
M.raw %>%
      ggplot(aes(x=year, y=M, group=as.character(Source))) +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(aes(colour=as.character(Source))) +
      labs(colour = "Keyrun") +
      guides(colour = guide_legend(nrow = 1)) +
      facet_wrap(~age,scales = 'free_y')
dev.off()

taf.png(file.path(paste0(prefix,"_M raw2_subset")))
# M by WR (facet) and Key run (colour)
subset(M.raw,Source %in% c(2016,2019,2023)) %>%
  ggplot(aes(x=year, y=M, group=as.character(Source))) +
  theme_bw() +
  theme(legend.position="bottom") +
  geom_line(aes(colour=as.character(Source))) +
  labs(colour = "Keyrun") +
  guides(colour = guide_legend(nrow = 1)) +
  facet_wrap(~age,scales = 'free_y')
dev.off()

taf.png(file.path(paste0(prefix,"_M decades")))
# M summarized by decade by WR (facet) and Key run (colour)
M.raw %>%
      mutate(decade = 10*floor(year/10)) %>% 
      group_by(Source, decade, age) %>% 
      summarise(M = mean(M, na.rm=TRUE)) %>% 
      group_by(Source, age) %>% 
      mutate(xend = lead(decade)) %>% 
      # View()
      ggplot(aes(x=decade, y=M, group=as.character(Source))) +
      theme_bw() +
      theme(legend.position="bottom") +
      # geom_line(aes(colour=as.character(Source))) +
      geom_segment(aes(xend=xend, yend=M, colour=as.character(Source))) +
      # scale_colour_manual(values=myPlotColors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(colour = "Keyrun") +
      expand_limits(y=0) +
      facet_wrap(~age,scales='free')
dev.off()

taf.png(file.path(paste0(prefix,"_M decades2")))
# M for ages 2-6 WR only, summarized by decade (x-axis) and Key run (colour)
M.raw %>%
      filter(year >= 1970) %>% 
      filter(age <= 6) %>% 
      mutate(decade = 10*floor(year/10)) %>% 
      mutate(WRgroup = ifelse(age %in% 0:1, "0-1","2-6")) %>% 
      group_by(Source, decade, WRgroup) %>% 
      summarise(M = mean(M, na.rm=TRUE)) %>% 
      group_by(Source, WRgroup) %>% 
      mutate(xend = lead(decade)) %>% 
      mutate(xend = ifelse(is.na(xend), Source, xend)) %>% 
      mutate(Mend = lead(M)) %>% 
      # View()
      ggplot(aes(x=decade, y=M, group=as.character(Source))) +
      theme_bw() +
      theme(legend.position="bottom") +
      # geom_line(aes(colour=as.character(Source))) +
      geom_segment(aes(xend=xend, yend=M, colour=as.character(Source))) +
      geom_segment(aes(x=xend, xend=xend, y=M, yend=Mend, colour=as.character(Source)), linetype = "dashed") +
      # scale_colour_manual(values=myPlotColors) +
      geom_point(aes(colour=as.character(Source))) +
      labs(colour = "Keyrun") +
      expand_limits(y=0) +
      facet_wrap(~WRgroup)
dev.off()

taf.png(file.path(paste0(prefix,"_M decades2_subset")))
# M for ages 2-6 WR only, summarized by decade (x-axis) and Key run (colour)
subset(M.raw,Source %in% c(2016,2019,2023)) %>%
  filter(year >= 1970) %>% 
  filter(age <= 6) %>% 
  mutate(decade = 10*floor(year/10)) %>% 
  mutate(WRgroup = ifelse(age %in% 0:1, "0-1","2-6")) %>% 
  group_by(Source, decade, WRgroup) %>% 
  summarise(M = mean(M, na.rm=TRUE)) %>% 
  group_by(Source, WRgroup) %>% 
  mutate(xend = lead(decade)) %>% 
  mutate(xend = ifelse(is.na(xend), Source, xend)) %>% 
  mutate(Mend = lead(M)) %>% 
  # View()
  ggplot(aes(x=decade, y=M, group=as.character(Source))) +
  theme_bw() +
  theme(legend.position="bottom") +
  # geom_line(aes(colour=as.character(Source))) +
  geom_segment(aes(xend=xend, yend=M, colour=as.character(Source))) +
  geom_segment(aes(x=xend, xend=xend, y=M, yend=Mend, colour=as.character(Source)), linetype = "dashed") +
  # scale_colour_manual(values=myPlotColors) +
  geom_point(aes(colour=as.character(Source))) +
  labs(colour = "Keyrun") +
  expand_limits(y=0) +
  facet_wrap(~WRgroup)
dev.off()

setwd('../..')

### ============================================================================
### load baseline
### ============================================================================
load(file.path(model.save,paste0('NSAS_HAWG2024_sf','.RData')))
NSH.sam.baseline  <- NSH.sam
NSH.baseline      <- NSH

load(file.path(model.save,paste0('NSAS_HAWG2024_sf_retro','.RData')))
# build mohn rho df
rho.retro.ssb <- mohns.rho(NSH.retro,ref.year = 2023,span = 10,type='ssb')
rho.retro.ssb$run <- 'retro_baseline'
rho.retro.ssb$type <- 'ssb'
rho.retro.fbar <- mohns.rho(NSH.retro,ref.year = 2023,span = 10,type='fbar')
rho.retro.fbar$run <- 'retro_baseline'
rho.retro.fbar$type <- 'fbar'
rho.retro.rec <- mohns.rho(NSH.retro,ref.year = 2023,span = 10,type='rec')
rho.retro.rec$run <- 'retro_baseline'
rho.retro.rec$type <- 'rec'

### ============================================================================
### profiling with 2023 M
### ============================================================================
prefix <- 'sf_scanM'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

load(file.path('../../model','NSAS_HAWG2024_sf_scanM.RData'))

res.scanM <- convert_scanM_results(res,NSH.sam.baseline,1974:2010,20)

df.nlogl    <- res.scanM$df.nlogl %>% mutate(min = (min(nlogl.norm) == nlogl.norm))
df.nlogl    <- df.nlogl %>% pivot_longer(!addM & !Mbar & !min,names_to ='type',values_to='value')

taf.png(file.path(paste0(prefix,"_profiling")))
ggplot(data=subset(df.nlogl,type == 'nlogl' | type == 'AIC'),aes(x=Mbar,y=value,col=as.factor(type)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'nlogl' | type == 'AIC') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  facet_wrap(~type,scales='free',ncol=1)+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2023)])))
dev.off()

taf.png(file.path(paste0(prefix,"_profiling addM")))
ggplot(data=subset(df.nlogl,type == 'nlogl' | type == 'AIC'),aes(x=addM,y=value,col=as.factor(type)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'nlogl' | type == 'AIC') & min==TRUE),aes(x=addM,y=value),color='red')+
  facet_wrap(~type,scales='free',ncol=1)+
  geom_vline(xintercept = 0)
dev.off()

setwd('../..')

### ============================================================================
### compare profiling between SMS runs
### ============================================================================
prefix <- 'sf_scanM_smsAll'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

runName <- paste0('NSAS_HAWG2024_',prefix)

flagFirst <- TRUE
for(SMSkeyRuns in SMSkeyRunsMat){
  load(file.path('../../model/',
                 paste0(runName,'_',SMSkeyRuns,"_sf.RData")))
  
  res.scanM <- convert_scanM_results(res,NSH.sam.baseline,1974:2010,20)
  df.nlogl.temp <- res.scanM$df.nlogl
  df.nlogl.temp$SMSkeyRuns  <- SMSkeyRuns
  
  if(flagFirst){
    df.nlogl            <- df.nlogl.temp
    flagFirst <- FALSE
  }else{
    df.nlogl                  <- rbind(df.nlogl,df.nlogl.temp)
  }
}

df.nlogl <- df.nlogl %>%
  group_by(SMSkeyRuns) %>%
  mutate(min = (min(nlogl.norm) == nlogl.norm)) %>%
  pivot_longer(!addM & !Mbar & !min & !SMSkeyRuns,names_to ='type',values_to='value')

saveTab <- df.nlogl[df.nlogl$min == TRUE,] %>% pivot_wider(names_from = type, values_from = value)
saveTab <- saveTab %>% select(addM,Mbar,SMSkeyRuns,nlogl,AIC,q,ssbAbs)
write.csv(x = saveTab,
          file = file.path(paste0('table_',runName,'.csv')),row.names = FALSE)

plotTab <- saveTab %>% pivot_longer(!addM & !SMSkeyRuns,names_to = 'type', values_to = 'value')
plotTab$run <- runName

taf.png(file.path(paste0(prefix,"_profiling")))
ggplot(data=subset(df.nlogl,type == 'nlogl' | type == 'AIC'),aes(x=Mbar,y=value,col=as.factor(SMSkeyRuns)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'nlogl' | type == 'AIC') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  facet_wrap(~type,scales='free',ncol=1)+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

taf.png(file.path(paste0(prefix,"_ssb ratio")))
ggplot(data=subset(df.nlogl,type == 'ssbRatio'),aes(x=Mbar,y=value,col=as.factor(SMSkeyRuns)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'ssbRatio') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  facet_wrap(~type,scales='free')+
  ylim(-0.5,0.5)+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

taf.png(file.path(paste0(prefix,"_q ratio")))
ggplot(data=subset(df.nlogl,type == 'q'),aes(x=Mbar,y=value,col=as.factor(SMSkeyRuns)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'q') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

taf.png(file.path(paste0(prefix,"_compare runs")))
p1 <- ggplot(data=subset(plotTab,type=='ssbAbs'),aes(x=as.factor(SMSkeyRuns),y=value))+
  geom_point()+
  ylim(1,2.5)+
  ylab('SSB')+
  ggtitle('SSB 20y span')+
  theme(axis.text.x = element_text(angle = 90))

p2 <- ggplot(data=subset(plotTab,type=='Mbar'),aes(x=as.factor(SMSkeyRuns),y=value))+
  geom_point()+
  ylim(0.2,0.5)+
  ylab('Mbar')+
  ggtitle('Mbar')+
  theme(axis.text.x = element_text(angle = 90))

p3 <- ggplot(data=subset(plotTab,type=='nlogl'),aes(x=as.factor(SMSkeyRuns),y=value))+
  geom_point()+
  ylab('nlogl')+
  ggtitle('neg log likelihood')+
  theme(axis.text.x = element_text(angle = 90))

p4 <- ggplot(data=subset(plotTab,type=='q'),aes(x=as.factor(SMSkeyRuns),y=value))+
  geom_point()+
  ylab('q')+
  ylim(0.8,2)+
  ggtitle('q HERAS(3-8)')+
  theme(axis.text.x = element_text(angle = 90))

grid.arrange(p1,p2,p3,p4)
dev.off()

setwd('../..')

### ============================================================================
### compare profiling between peels
### ============================================================================
prefix <- 'sf_scanM_peels'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

runName <- paste0('NSAS_HAWG2024_',prefix)

flagFirst <- TRUE
nPeels <- 10
endYear <- 2023

for(iPeel in 1:nPeels){
  ### ============================================================================
  ### load results
  ### ============================================================================
  
  load(file.path('../../model',paste0(runName,'_peel',endYear-(iPeel-1),"_sf.RData")))
  
  ### ============================================================================
  ### formatting
  ### ============================================================================
  
  res.scanM <- convert_scanM_results(res,NSH.sam.baseline,1974:2010,20)
  df.nlogl.temp <- res.scanM$df.nlogl
  df.nlogl.temp$peel  <- endYear-(iPeel-1)
  #df.nlogl.temp$mType     <- iType
  
  if(flagFirst){
    df.nlogl            <- df.nlogl.temp
    flagFirst <- FALSE
  }else{
    df.nlogl                  <- rbind(df.nlogl,df.nlogl.temp)
  }
}

df.nlogl <- df.nlogl %>%
  group_by(peel) %>%
  mutate(min = (min(nlogl.norm) == nlogl.norm)) %>%
  pivot_longer(!addM & !Mbar & !min & !peel,names_to ='type',values_to='value')

saveTab <- df.nlogl[df.nlogl$min == TRUE,] %>% pivot_wider(names_from = type, values_from = value)
saveTab <- saveTab %>% select(addM,Mbar,peel,nlogl,AIC,q,ssbAbs)
write.csv(x = saveTab,
          file = file.path(paste0('table_',runName,'.csv')),row.names = FALSE)

plotTab <- saveTab %>% pivot_longer(!addM & !peel,names_to = 'type', values_to = 'value')
plotTab$run <- runName


# plot metrics by peel
taf.png(file.path(paste0(prefix,"_metrics")))
ggplot(data=plotTab,aes(x=as.factor(peel),y=value))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~type,scales = 'free')
dev.off()

# normalized log likelihood and AIC
taf.png(file.path(paste0(prefix,"_profiling")))
ggplot(data=subset(df.nlogl,type == 'nlogl' | type == 'AIC'),aes(x=Mbar,y=value,col=as.factor(peel)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'nlogl' | type == 'AIC') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  facet_wrap(~type,scales='free',ncol=1)+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

# SSB ratio
taf.png(file.path(paste0(prefix,"_ssb ratio")))
ggplot(data=subset(df.nlogl,type == 'ssbRatio'),aes(x=Mbar,y=value,col=as.factor(peel)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'ssbRatio') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  facet_wrap(~type,scales='free')+
  ylim(-0.5,0.5)+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

# q ratio
taf.png(file.path(paste0(prefix,"_q ratio")))
ggplot(data=subset(df.nlogl,type == 'q'),aes(x=Mbar,y=value,col=as.factor(peel)))+
  geom_line()+
  geom_point(data = subset(df.nlogl,(type == 'q') & min==TRUE),aes(x=Mbar,y=value),color='red')+
  geom_vline(xintercept = yearMeans(quantMeans(NSH.baseline@m[,ac(1974:2010)])))
dev.off()

# retro profile
load(file.path('../../model','NSAS_HAWG2024_sf_scanM_peels.RData'))
rho.prof.ssb <- mohns.rho(NSH.sams,ref.year = 2023,span = 10,type='ssb')
rho.prof.ssb$run <- runName
rho.prof.ssb$type <- 'ssb'
rho.prof.fbar <- mohns.rho(NSH.sams,ref.year = 2023,span = 10,type='fbar')
rho.prof.fbar$run <- runName
rho.prof.fbar$type <- 'fbar'
rho.prof.rec <- mohns.rho(NSH.sams,ref.year = 2023,span = 10,type='rec')
rho.prof.rec$run <- runName
rho.prof.rec$type <- 'rec'

df.rho <- rbind(rho.prof.ssb,
                rho.prof.fbar,
                rho.prof.rec,
                rho.retro.ssb,
                rho.retro.fbar,
                rho.retro.rec)

taf.png(file.path(paste0(prefix,"_retro M")))
print(plot(NSH.sams))
dev.off()

# mohn rho
taf.png(file.path(paste0(prefix,"_mohn rho")))
ggplot(data=df.rho,aes(x=year,y=rho,col=run))+
  geom_line()+
  facet_wrap(~type,ncol = 1,scales = 'free')
dev.off()

# plot metrics by scenario
taf.png(file.path(paste0(prefix,"_mohn rho")))
p2 <- ggplot(data=subset(plotTab,type=='Mbar'),aes(x=as.factor(peel),y=value))+
  geom_point()+
  ylim(0.2,0.5)+
  ylab('Mbar')+
  ggtitle('Mbar')+
  theme(axis.text.x = element_text(angle = 90))

p3 <- ggplot(data=subset(plotTab,type=='AIC'),aes(x=as.factor(peel),y=value))+
  geom_point()+
  ylab('nlogl')+
  ggtitle('neg log likelihood')+
  theme(axis.text.x = element_text(angle = 90))

p4 <- ggplot(data=subset(plotTab,type=='q'),aes(x=as.factor(peel),y=value))+
  geom_point()+
  ylab('q')+
  ylim(0.8,2)+
  ggtitle('q HERAS(3-8)')+
  theme(axis.text.x = element_text(angle = 90))

grid.arrange(grobs = list(p2,p3,p4),
                  layout_matrix = rbind(c(NA,1),
                                        c(2, 3)))
dev.off()
