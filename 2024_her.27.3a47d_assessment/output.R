## Extract results of interest, write TAF output tables

## Before: results_sf.RData (model)
## After:  fatage.csv, natage.csv, results_sf.RData, summary.csv (output)

rm(list=ls())

library(icesTAF)
library(tidyverse)

taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)

runName <- 'NSAS_HAWG2024'

mkdir("output")

data.save       <- file.path(".", "data")
model.save      <- file.path(".",'model')
output.save     <- file.path(".", "output")

terminal.year <- 2023

#cp("model/results_sf.RData", "output")
#cp("model/results_mf.RData", "output")

load(file.path(model.save,paste0(runName,'_sf.RData')))
load(file.path(model.save,paste0(runName,'_sf_retro.RData')))
load(file.path(model.save,paste0(runName,'_mf.RData')))
load(file.path(model.save,paste0(runName,'_mf_retro.RData')))


flagFirst <- TRUE
for(idx in 1:length(NSH.retro)){
  q <- catchabilities(NSH.retro[idx])
  
  if(flagFirst){
    q.all <- q
    flagFirst <- FALSE
  }else{
    q.all <- rbind(q.all,q)
  }
}

# windows()
# ggplot(subset(q.all,fleet == 'HERAS'),aes(x=age,y=value,col=name))+
#   geom_line()
# 
# windows()
# ggplot(subset(q.all,fleet == 'IBTS-Q1'),aes(x=age,y=value,col=name))+
#   geom_point()
# 
# windows()
# ggplot(subset(q.all,fleet == 'IBTS-Q3'),aes(x=age,y=value,col=name))+
#   geom_line()

### ============================================================================
### Single fleet outputs
### ============================================================================
write.taf(flr2taf(NSH.sam@stock.n), "output/natage_sf.csv")
write.taf(flr2taf(NSH.sam@harvest), "output/fatage_sf.csv")

components.df <- as.data.frame(components(NSH.sam))
components.df <- components.df[,1:dim(components.df)[2]]
colnames(components.df) <- ac(1972:range(NSH)['maxyear'])

write.taf(t(components.df), "output/components.csv")
write.taf(as.data.frame(NSH.tun$IBTS0@index), "output/IBTS0.csv")

IBTSQ1      <- read.table(file.path(".", "bootstrap",'data','DATRAS_IBTSQ1.csv'), sep=",", header = TRUE) # read raw indices instead of standardized ones
IBTSQ1      <- IBTSQ1[IBTSQ1$IndexArea == "NS_Her",]
IBTSQ1 <- data.frame(cbind(IBTSQ1$Year, IBTSQ1$Age_1))
names(IBTSQ1) <- c("year","data")
write.taf(IBTSQ1, "output/IBTSQ1.csv")

NSH.tun$IBTS0

## Summary
vlu <- c("value", "lbnd", "ubnd")

summary <- data.frame(
  rec(NSH.sam)$year,
  rec(NSH.sam)[vlu],
  tsb(NSH.sam)[vlu],
  ssb(NSH.sam)[vlu],
  catch(NSH.sam)[vlu],
  fbar(NSH.sam)[vlu],
  c(catch(NSH), NA),
  c(sop(NSH), NA), row.names=NULL)

names(summary) <- c(
  "Year",
  "Rec",   "Rec_lo",   "Rec_hi",
  "TSB",   "TSB_lo",   "TSB_hi",
  "SSB",   "SSB_lo",   "SSB_hi",
  "Catch", "Catch_lo", "Catch_hi",
  "Fbar",  "Fbar_lo",  "Fbar_hi",
  "Landings", "SOP")
write.taf(summary, "output/summary_sf.csv")

df.mr <- data.frame(matrix(ncol = 3, nrow = 0))


df.mr <- rbind(df.mr,cbind(seq((terminal.year-10),terminal.year),
                           t(t(mohns.rho(NSH.retro,span=10,ref.year=terminal.year,type="ssb")$rho)),
                           rep('ssb',10,1)))
df.mr <- rbind(df.mr,cbind(seq((terminal.year-10),terminal.year),
                           t(t(mohns.rho(NSH.retro,span=10,ref.year=terminal.year,type="fbar")$rho)),rep('fbar',10,1)))
df.mr <- rbind(df.mr,cbind(seq((terminal.year-10),terminal.year),
                           t(t(mohns.rho(NSH.retro,span=10,ref.year=terminal.year,type="rec")$rho)),rep('rec',10,1)))

df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(NSH.retro,span=5,ref.year=terminal.year,type="ssb")$rho),'ssb'))
df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(NSH.retro,span=5,ref.year=terminal.year,type="fbar")$rho),'fbar'))
df.mr <- rbind(df.mr,cbind('av_5y',mean(mohns.rho(NSH.retro,span=5,ref.year=terminal.year,type="rec")$rho),'rec'))

colnames(df.mr) <- c('year','mohn_rho','var')

df.mr <- df.mr %>% pivot_wider(names_from = var, values_from = mohn_rho)

write.taf(df.mr, "output/mohn_rho_sf.csv")

### ============================================================================
### Multi fleet outputs
### ============================================================================

## N at age, F at age
write.taf(flr2taf(NSH3f.sam@stock.n), "output/natage_mf.csv")

write.taf(flr2taf(NSH3f.sam@harvest[,,,,'A']), "output/fatage_A_mf.csv")
write.taf(flr2taf(NSH3f.sam@harvest[,,,,'BD']), "output/fatage_BD_mf.csv")
write.taf(flr2taf(NSH3f.sam@harvest[,,,,'C']), "output/fatage_C_mf.csv")

# Function flr2taf does not handle multifleet yet! MP, 23/3/2021
# write.taf(flr2taf(NSH3f.sam@harvest), "output/fatage_mf.csv")


# summary of multifleet
summary <- data.frame(
  rec(NSH3f.sam)$year,
  rec(NSH3f.sam)[vlu],
  tsb(NSH3f.sam)[vlu],
  ssb(NSH3f.sam)[vlu],
  catch(NSH3f.sam)[vlu],
  fbar(NSH3f.sam)[vlu],
  c(areaSums(catch(NSH3f)),NA),
  # c(sop(NSHs3$residual), NA), 
  row.names=NULL)

names(summary) <- c(
  "Year",
  "Rec",   "Rec_lo",   "Rec_hi",
  "TSB",   "TSB_lo",   "TSB_hi",
  "SSB",   "SSB_lo",   "SSB_hi",
  "Catch", "Catch_lo", "Catch_hi",
  "Fbar",  "Fbar_lo",  "Fbar_hi",
  "Landings")
write.taf(summary, "output/summary_mf.csv")

