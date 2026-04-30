# ==============================================================================
# Reference points estimation
# North Sea herring at interbenchmark meeting 2021
# Sensititivity analysis with ricker and segreg
#
# 22/02/2018 Martin Pastoors
# 17/03/2018 Checked and updated folders after cleanup of github
# 07/06/2021 Prepared for IBP 2021 meeting
# 08/06/2021 Adapted to new reference points systematics (David Miller)
# 18/08/2021 Final after IBP 2021 meeting
# 09/05/2022 Testing with 2022 assessment

# Steps:
# 1. Get estimate of Blim at breakpoint using the whole time series and calculate Bpa from it
# 2. parameterize the segreg model with Blim breakpoint and (roughly) geomean rec above this
# 3. truncate the NSH object to start after 2002
# 4. fit the stock recruitment model(s)
# 5. Get Flim and thereby Fpa. Run EqSim with no MSY Btrigger (i.e. run EqSim with Btrigger=0), and Fcv=Fphi=0
# 6. Run EqSim with assessment error but no MSY Btrigger (i.e. run EqSim with Btrigger=0),
#    to get initial FMSY ; if this initial FMSY value is > Fpa, reduce it to Fpa
# 7. Check if FMSY is precautionary, so do a scan on Fp05. If Fmsy is larger than Fp05, reduce to Fp05
# 8. final set of reference points
#
# ==============================================================================

#rm(list=ls())

library(icesTAF)

taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)

# follow the order of loading (1) FLCore and (2) msy
# as SR functions have same name but different formulation
taf.library(msy)

# taf.session()

library(methods)
library(tidyverse)

mkdir("refpoints")
mkdir(file.path("refpoints","input"))
mkdir(file.path("refpoints","output"))
mkdir(file.path("refpoints","plots"))

# source("Refpoints functions.R")

# settings

Fcv       = 0.16            #From IBPNSAS 2021; 12 years
Fphi      = 0.47            #From IBPNSAS 2021; 12 years

basename = gsub("\\.RData","",filename)

basename_main = gsub("\\_retro.RData","",filename)

# load assessment data
load(file.path("model", filename))

load(file.path("model", paste0(basename_main,'.RData')))
# load(file.path("refpoints","input", "WKPELA2018_NSH_final.RData"))

# temporary fix
# NSH.sam <- NSH.samBind1
NSHbase <- NSH
refpts  <- list()
for(iRetro in names(NSH.retro)){
  NSH <- NSHbase + NSH.retro[[iRetro]]
  # get max year and set ranges
  maxyear <- as.integer(NSH@range["maxyear"])
  bio.years = c((maxyear-9),maxyear)
  sel.years = c((maxyear-9),maxyear)


  # --------------------------------------------------------------------------------------------------------------  
  # 1. Get estimate of Blim at breakpoint using the whole time series and calculate Bpa from it
  # --------------------------------------------------------------------------------------------------------------  
  FIT_segregBlim <- eqsr_fit(NSH,nsamp=2000, models = "Segreg", rshift=1, remove.years = 1979:1990)

  eqsr_plot(FIT_segregBlim, n=2e4, ggPlot=TRUE)
  # ggsave(filename=file.path("refpoints","plots",paste0(basename,"_eqsr_segregblim.jpg")), device="jpeg")

  Blim <- round(FIT_segregBlim$sr.det$b/1e5)*1e5  # 796 kT = 800 kT
  Blim <- round(FIT_segregBlim$sr.det$b)  

  #Now calculate the uncertainty in SSB in terminal year. We need the sd that belongs to log(ssb) to calculate Bpa
  logssb    <- subset(ssb(NSH.sam),year==maxyear)
  sdmin     <- function(sdestim){
    return(abs(0.025 - dnorm(log(logssb$lbnd),log(logssb$value),sdestim)))}
  sdSSB     <- optimize(sdmin,interval=c(1e-4,0.2))$minimum
  # sdSSB     <- max(sdSSB, 0.2)

  Bpa       <- round(Blim * exp(1.645*sdSSB))  


  # --------------------------------------------------------------------------------------------------------------  
  # 2. parameterize the segreg model with Blim breakpoint and (roughly) geomean rec above this
  # --------------------------------------------------------------------------------------------------------------  
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, 
                                              ab$a * Blim, 
                                              ab$a * ssb))

  # --------------------------------------------------------------------------------------------------------------  
  # 3. truncate the NSH object
  # --------------------------------------------------------------------------------------------------------------  
  NSHtrunc <- trim(NSH, year=2002:maxyear)
  # NSHtrunc <- trim(NSH, year=2002:2016)

  # --------------------------------------------------------------------------------------------------------------  
  # 4. fit the stock recruitment model(s)
  # --------------------------------------------------------------------------------------------------------------  
  # FIT <- eqsr_fit(NSHtrunc, nsamp = 2000, models = c("Ricker", "SegregBlim"), rshift=1)
  FIT <- eqsr_fit(NSHtrunc, nsamp = 2000, models = c("SegregBlim"), rshift=1)
  # FIT <- eqsr_fit(NSHtrunc, nsamp = 2000, models = c("Segreg"), rshift=1)
  # FIT <- eqsr_fit(NSHtrunc, nsamp = 2000, models = c("Cadigan"), rshift=1)
  # FIT <- eqsr_fit(NSHtrunc, nsamp = 100, models = c("Ricker", "SegregBlim"), rshift=1)

  eqsr_plot(FIT,n=2e4, ggPlot=TRUE)
  # ggsave(filename=file.path("refpoints","plots",paste0(basename,"_eqsr_fit_segregblim.jpg")), device="jpeg")

  # round(FIT$sr.det$b) / 1.4



  # --------------------------------------------------------------------------------------------------------------  
  # 5. Get Flim and thereby Fpa. Run EqSim with no MSY Btrigger (i.e. run EqSim with Btrigger=0), and Fcv=Fphi=0
  # --------------------------------------------------------------------------------------------------------------  
  SIM1 <- eqsim_run(FIT,
                   bio.years        = bio.years,
                   bio.const        = FALSE,
                   sel.years        = sel.years,
                   sel.const        = FALSE,
                   recruitment.trim = c(3, -3),
                   Fcv              = 0,
                   Fphi             = 0,
                   Blim             = Blim,
                   Bpa              = Bpa,
                   Btrigger         = 0,
                   Fscan            = seq(0,0.80,len=40),
                   # Fscan            = c(seq(0,0.2,by=0.001), seq(0.21,0.8,by=0.01)),
                   verbose          = TRUE,
                   extreme.trim     = c(0.01,0.99))

  # msy::eqsim_plot(SIM1, ggPlot=TRUE)

  Flim      <- SIM1$Refs2["catF","F50"]   # MP: 0.341

  # Now calculate the uncertainty in F in terminal year. We need the sd that belongs to log(F) to calculate Fpa
  # logfbar   <- subset(fbar(NSH.sam),year==maxyear)
  # sdmin     <- function(sdestim){
  #              return(abs(0.025 - dnorm(log(logfbar$lbnd),log(logfbar$value),sdestim)))}
  # sdF       <- optimize(sdmin,interval=c(1e-4,0.2))$minimum

  #DM: old method commented out
  #Fpa       <- Flim * exp(-1.645*sdF) #0.294  MP: 0.298     

  # --------------------------------------------------------------------------------------------------------------  
  # 6. Run EqSim with assessment error but no MSY Btrigger (i.e. run EqSim with Btrigger=0),
  #    to get initial FMSY ; 
  #   (#DM: not yet, check Fp05=Fpa later) if this initial FMSY value is > Fpa, reduce it to Fpa
  # --------------------------------------------------------------------------------------------------------------  
  SIM2 <- eqsim_run(FIT,
                    bio.years = bio.years,
                    bio.const = FALSE,
                    sel.years = sel.years,
                    sel.const = FALSE,
                    recruitment.trim = c(3, -3),
                    Fcv       = Fcv,            
                    Fphi      = Fphi,            
                    Blim      = Blim,
                    Bpa       = Bpa,
                    Btrigger  = 0,
                    Fscan     = seq(0,0.80,len=40),
                    # Fscan  = c(seq(0,0.2,by=0.001), seq(0.21,0.8,by=0.01)),
                    verbose   = TRUE,
                    extreme.trim=c(0.01,0.99))

  Fmsy      <- SIM2$Refs2["lanF","medianMSY"] #0.275
  #DM: not yet, check Fp05=Fpa later
  #Fmsy      <- ifelse(Fmsy>Fpa,Fpa,Fmsy)
  # msy::eqsim_plot(SIM2)



  # Select MSY Btrigger   (from schematic guidelines: yes, yes, no -> 5th percentile of MSYBtrigger
  MSYBtrigger <- SIM2$Refs2["catB","F05"]  # MP 1396 kT
  MSYBtrigger <- round(MSYBtrigger) # rounding

  # --------------------------------------------------------------------------------------------------------------  
  # 7. Check if FMSY is precautionary, so do a scan on Fp05. If Fmsy is larger than Fp05, reduce to Fp05
  # --------------------------------------------------------------------------------------------------------------  
  SIM3 <- eqsim_run(FIT,
                    bio.years = bio.years,
                    bio.const = FALSE,
                    sel.years = sel.years,
                    sel.const = FALSE,
                    recruitment.trim = c(3, -3),
                    Fcv       = Fcv,
                    Fphi      = Fphi,
                    Blim      = Blim,
                    Bpa       = Bpa,
                    Btrigger  = MSYBtrigger,
                    Fscan     = seq(0,0.80,len=40),
                    # Fscan     = seq(0.18,0.40,0.01),
                    # Fscan     = seq(0.18,0.5,0.01),
                    # Fscan     = c(seq(0,0.2,by=0.001), seq(0.21,0.8,by=0.01)),
                    verbose   = TRUE,
                    extreme.trim=c(0.01,0.99))

  # msy::eqsim_plot(SIM3, ggPlot=TRUE)

  # If the precautionary criterion (FMSY < Fp.05) evaluated is not met, then FMSY should be reduced to  Fp.05. 
  Fp05      <- SIM3$Refs2["catF","F05"] # MP: 0.256
  #DM: define new Fpa here
  Fpa <- Fp05
  #DM: if Fpa > Flim, then Flim will be undefined
  if (Fpa>Flim) Flim <- NA

  propFmsy  <- subset(SIM3$pProfile, round(Ftarget, 2)==round(Fmsy,2) & variable=="Blim")$value
  if (Fmsy > Fp05) {Fmsy <- Fp05}

  # --------------------------------------------------------------------------------------------------------------  
  # 8. final set of reference points
  # --------------------------------------------------------------------------------------------------------------  
  #DM - use ICES rounding
  if (Flim<0.2) Flim   <- round(Flim,3) else Flim   <- round(Flim,2)
  if (Fpa<0.2) Fpa     <- round(Fpa,3) else Fpa   <- round(Fpa,2)
  if (Fmsy<0.2) Fmsy   <- round(Fmsy,3) else Fmsy   <- round(Fmsy,2)

  #Flim   <- round(Flim,2)
  #Fpa    <- round(Fpa,2)
  #Fmsy   <- round(Fmsy, 2)
  refpts[[iRetro]] <- data.frame(Flim       = Flim,
                       Fpa        = Fpa,
                       Fmsy       = Fmsy,
                       Fp05       = Fp05,
                       Blim       = Blim,
                       Bpa        = Bpa,
                       MSYBtrigger= MSYBtrigger,
                       Fcv        = Fcv,
                       Fphi       = Fphi)
}

save.image(file = paste0(basename,'_refpoints.RData'))

refpts <- do.call(rbind,
                  lapply(1:length(refpts),
                         function(x){cbind(yearrun=c(2023:2013)[x],
                                           refpts[[x]])}))

write.taf(refpts,file = file.path("refpoints","output",paste0(basename,'.csv')))


# print(refpts)
pander::pandoc.table(refpts, 
                     style        = "simple",
                     split.tables = 200, 
                     split.cells  = c(rep(7,10)),
                     justify      = "right",
                     missing      =" ",
                     big.mark     = '', 
                     round        = c(2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0))

save(NSH, NSHtrunc, FIT_segregBlim, FIT, 
     SIM1, SIM2, SIM3,
     refpts, 
     file=file.path("refpoints","output", paste0(basename,"_refpoints_segregblim.RData")))


# load(file=file.path("refpoints","output", paste0(basename,"_refpoints_segregblim.RData")))

# Temporary comparison

# IBP <- cbind(rec(NSH.sam)$year,
#              rec(NSH.sam)$value, 
#              ssb(NSH.sam)$value,  
#              fbar(NSH.sam)$value,
#              c(catch(NSH),NA)) 
# colnames(IBP) <- c("year","rec","ssb","fbar","catch")
# IBP <- as_tibble(IBP) %>% mutate(source="IBP")

# load(file.path("model","assessment","WKPELA2018_final.RData"))
# WKPELA <- cbind(rec(NSH.sam)$year,
#              rec(NSH.sam)$value, 
#              ssb(NSH.sam)$value,  
#              fbar(NSH.sam)$value,
#              c(catch(NSH),NA)) 
# colnames(WKPELA) <- c("year","rec","ssb","fbar","catch")
# WKPELA <- as_tibble(WKPELA) %>% mutate(source="WKPELA")

# bind_rows(IBP, WKPELA) %>% 
#   pivot_longer(names_to="variable", values_to="data", rec:catch) %>% 
#   filter(year >= 2000) %>% 
#   
#   ggplot(aes(x=year, y=data)) +
#   theme_bw() +
#   geom_line(aes(colour=source), size=1) +
#   facet_wrap(~variable, scales = "free_y")

# sessionInfo()
# 
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_Netherlands.1252  LC_CTYPE=English_Netherlands.1252    LC_MONETARY=English_Netherlands.1252
# [4] LC_NUMERIC=C                         LC_TIME=English_Netherlands.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] forcats_0.5.1         stringr_1.4.0         dplyr_1.0.6           purrr_0.3.4           readr_1.4.0          
# [6] tidyr_1.1.3           tibble_3.1.2          ggplot2_3.3.3         tidyverse_1.3.0       msy_0.1.19           
# [11] FLSAM_2.1.1           stockassessment_0.5.4 FLCore_2.6.15.9023    iterators_1.0.13      lattice_0.20-41      
# [16] icesTAF_3.6.0        
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.2         jsonlite_1.7.2     splines_4.0.3      ellipse_0.4.2      modelr_0.1.8       assertthat_0.2.1  
# [7] stats4_4.0.3       pander_0.6.3       cellranger_1.1.0   yaml_2.2.1         ggrepel_0.9.1      pillar_1.6.1      
# [13] backports_1.2.1    glue_1.4.2         digest_0.6.27      RColorBrewer_1.1-2 rvest_0.3.6        colorspace_2.0-0  
# [19] htmltools_0.5.1.1  Matrix_1.2-18      plyr_1.8.6         pkgconfig_2.0.3    broom_0.7.6        haven_2.3.1       
# [25] scales_1.1.1       mgcv_1.8-33        generics_0.1.0     farver_2.0.3       ellipsis_0.3.2     withr_2.4.1       
# [31] TMB_1.7.20         cli_2.5.0          magrittr_2.0.1     crayon_1.4.1       readxl_1.3.1       evaluate_0.14     
# [37] fs_1.5.0           fansi_0.4.2        nlme_3.1-149       MASS_7.3-53        xml2_1.3.2         tools_4.0.3       
# [43] hms_1.0.0          lifecycle_1.0.0    munsell_0.5.0      reprex_1.0.0       compiler_4.0.3     rlang_0.4.10      
# [49] grid_4.0.3         rstudioapi_0.13    labeling_0.4.2     rmarkdown_2.6      gtable_0.3.0       DBI_1.1.1         
# [55] roxygen2_7.1.1     reshape2_1.4.4     R6_2.5.0           lubridate_1.7.9.2  knitr_1.31         utf8_1.1.4        
# [61] stringi_1.5.3      parallel_4.0.3     Rcpp_1.0.6         vctrs_0.3.8        dbplyr_2.1.0       tidyselect_1.1.0  
# [67] xfun_0.20
