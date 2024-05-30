scanM_fun <- function(  NSH,
                        NSH.tun,
                        NSH.ctrl,
                        mOrig,
                        iM){#stk0.NSH.sam,
  # sequence in case we don't have a convergence
  stk0.addM <- seq(-0.1,0.6,0.1)
  
  #stk0.NSH.sam.init <- stk0.NSH.sam
  NSH.init          <- NSH
  print(iM)
  
  NSH@m <- mOrig + iM
  res   <- try(FLSAM(NSH,NSH.tun,NSH.ctrl))#,starting.values = stk0.NSH.sam.init
  
  countLoop <- 1
  while(class(res)=="try-error" & countLoop <= length(stk0.addM)){
    NSH.init@m <- mOrig + stk0.addM[countLoop]
    #stk0.NSH.sam   <- try(FLSAM(NSH.init,NSH.tun,NSH.ctrl,starting.values = stk0.NSH.sam.init))
    res   <- try(FLSAM(NSH,NSH.tun,NSH.ctrl))#,starting.values = stk0.NSH.sam
    
    countLoop <- countLoop + 1
  }
  
  if(countLoop > length(stk0.addM) & class(res)=="try-error"){
    res <- new("FLSAM")
  }
  
  return(list(res.sam=res,addM=iM,m=NSH@m))
}

run_scanM <- function(addM,NSH,NSH.tun,NSH.ctrl,mOrig){#stk0.NSH.sam,
  require(doParallel)
  ncores <- detectCores()-1
  ncores <- ifelse(length(addM)<ncores,length(addM),ncores)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  clusterEvalQ(cl,library(icesTAF))
  clusterEvalQ(cl,source('utilities_model.R'))
  clusterEvalQ(cl,taf.library(FLCore))
  clusterEvalQ(cl,taf.library(stockassessment))
  clusterEvalQ(cl,taf.library(FLSAM))
  res     <- foreach(iM = addM) %dopar% scanM_fun(NSH,NSH.tun,NSH.ctrl,mOrig,iM)#,stk0.NSH.sam
  stopCluster(cl)
  detach("package:doParallel",unload=TRUE)
  detach("package:foreach",unload=TRUE)
  
  NSH.sams  <- new("FLSAMs")
  NSH.m     <- new("FLQuants")
  for(i in 1:length(res)){
    NSH.sams[[ac(res[[i]]$addM)]] <- res[[i]]$res.sam
    NSH.m[[ac(res[[i]]$addM)]]    <- res[[i]]$m
  }
  
  return(list(NSH.sams=NSH.sams,NSH.m=NSH.m))
}


#' @title Converts FLStock object to rbya
#'
#' @description The rbya (results by year and age) is a long \code{data.frame}
#' popular in Reykjavik
#'
#' @export
#'
#' @param x An FLStock object
#' @param scale scaler used on abundance values (stock in numbers, catches etc)
#' @param project A boolean, if TRUE (default), propagates terminal stock numbers
#' forward by one year (into the assessment year). Note that the weights in the
#' assessment year are the same as in the terminal year.
#' @param plusgroup A boolean, if TRUE (default), last age group is a plus group.
#' Only used if project is TRUE.

flstock_to_rbya <- function(x, scale=1, project = TRUE, plusgroup = TRUE)
{
  
  y    <- reshape2::melt(x@stock.n@.Data,value.name = "n")[,c("year","age","n")]
  y$n  <- y$n/scale
  y$f  <- reshape2::melt(x@harvest@.Data)[,c("value")]
  # if(class(x) != "FLSAM") {  # This may be needed
  y$oC <- reshape2::melt(x@catch.n@.Data)[,c("value")]/scale
  y$cW <- reshape2::melt(x@catch.wt@.Data)[,c("value")]
  y$sW <- reshape2::melt(x@stock.wt@.Data)[,c("value")]
  y$oD  = reshape2::melt(x@discards.n@.Data)[,c("value")]/scale
  y$dW  = reshape2::melt(x@discards.wt@.Data)[,c("value")]
  y$oL  = reshape2::melt(x@landings.n@.Data)[,c("value")]/scale
  y$lW  = reshape2::melt(x@landings.wt@.Data)[,c("value")]
  y$mat = reshape2::melt(x@mat@.Data)[,c("value")]
  y$pF  = reshape2::melt(x@harvest.spwn@.Data)[,c("value")]
  y$pM  = reshape2::melt(x@m.spwn@.Data)[,c("value")]
  y$m   = reshape2::melt(x@m@.Data)[,c("value")]
  
  # propagate stock forward
  if (project) {
    y2 <- y[y$year == max(y$year),]
    y2$year <- y2$year + 1
    y2$n <- y2$n * exp(-(y2$m + y2$f))
    if(plusgroup) {
      y2$n[(nrow(y2)-1)] <- y2$n[(nrow(y2)-1)] + y2$n[nrow(y2)]
    }
    y2$n <- c(NA, y2$n[2:length(y2$n)])
    y2$f <- y2$oC <- y2$oD <- y2$oL <- NA
    
    y <- rbind(y, y2)
  }
  
  return(dplyr::as_tibble(y))
  
}

convert_scanM_results <- function(res,NSH.sam.baseline,yearRangeM,yearSpanSSB){
  addM <- as.numeric(names(res$NSH.m))
  
  Mbar              <- array(NA, dim=c(length(addM),1))
  nlogl             <- array(NA, dim=c(length(addM),1))
  AIC               <- array(NA, dim=c(length(addM),1))
  ssbRatio          <- array(NA, dim=c(length(addM),1))
  ssbAbs            <- array(NA, dim=c(length(addM),1))
  q                 <- array(NA, dim=c(length(addM),1))

  flagFirst <- TRUE
  for(idx in 1:length(addM)){
    Mbar[idx] <- yearMeans(quantMeans(res$NSH.m[[idx]][,ac(yearRangeM)]))
    
    # only loop on those that converged
    if(dims(res$NSH.sams[[idx]])$year != 1){
      # extract log likelihood
      nlogl[idx]  <- nlogl(res$NSH.sams[[idx]])
      AIC[idx]    <- AIC(res$NSH.sams[[idx]])
      
      # extract SSB relative to baseline
      ssb.baseline  <- ssb(NSH.sam.baseline)
      ssb.baseline  <- ssb.baseline[  ssb.baseline$year <= range(res$NSH.sams[[idx]])['maxyear'] & 
                                        ssb.baseline$year >= range(res$NSH.sams[[idx]])['minyear'],]
      
      ssb.current               <- ssb(res$NSH.sams[[idx]])
      ssb.current$baselineRatio <- (ssb.baseline$value-ssb.current$value)/ssb.baseline$value
      
      ssbRatio[idx] <- mean(ssb.current$baselineRatio)
      
      # extract SSB level
      temp        <- subset(ssb.current,year >= (range(res$NSH.sams[[idx]])['maxyear']-yearSpanSSB))
      ssbAbs[idx] <- mean(temp$value)*1e-6
      
      # extract catchability HERAS age 3+
      qMat <- catchabilities(res$NSH.sams[[idx]])
      qMat <- qMat[qMat$fleet == 'HERAS' & qMat$age>=3,]
      q[idx] <- mean(qMat$value)
      
      # extract parameters
      params.select <- unique(colnames(res$NSH.sams[[idx]]@vcov))
      
      param.table <- res$NSH.sams[[idx]]@params
      param.table <- param.table[which(param.table$name %in% params.select),]
      
      for(idxParam in params.select){
        nvar <- length(which(param.table$name == idxParam))
        param.table$name[param.table$name == idxParam] <- paste(param.table$name[param.table$name == idxParam],1:nvar,sep='.')
      }
      
      param.table     <- param.table %>% pivot_longer(!name,names_to ='type',values_to='value',names_repair = "unique")
      param.table$run <- addM[idx]

      if(flagFirst){
        param.table.full <- param.table
        flagFirst <- FALSE
      }else{
        param.table.full <- rbind(param.table.full,param.table)
      }
    }
  }
  
  # filter NA (no convergence)
  Mbar              <- Mbar[!is.na(nlogl)]
  addM              <- addM[!is.na(nlogl)]
  ssbRatio          <- ssbRatio[!is.na(nlogl)]
  ssbAbs            <- ssbAbs[!is.na(nlogl)]
  q                 <- q[!is.na(nlogl)]
  AIC               <- AIC[!is.na(nlogl)]
  nlogl             <- nlogl[!is.na(nlogl)]
  
  # filter outliers
  idxFilt <- which(nlogl<(quantile(nlogl, 0.5)+100) & nlogl>(quantile(nlogl, 0.5)-100))
  
  Mbar              <- Mbar[idxFilt]
  addM              <- addM[idxFilt]
  ssbRatio          <- ssbRatio[idxFilt]
  ssbAbs            <- ssbAbs[idxFilt]
  q                 <- q[idxFilt]
  nlogl             <- nlogl[idxFilt]
  AIC               <- AIC[idxFilt]

  nlogl.norm    <- nlogl-mean(nlogl)
  
  # convert to df
  df.nlogl  <- as.data.frame(cbind(addM,Mbar,nlogl,AIC,nlogl.norm,ssbRatio,q,ssbAbs))
  
  return(list(df.nlogl=df.nlogl,
              param.table = param.table.full))
}

comp_sam_params <- function(NSH.sams){
  
  runNames <- names(NSH.sams)
  
  flagFirst <- TRUE
  for(run in runNames){
    if(dims(NSH.sams[[run]])$year != 1){
      # stock metrics
      ssb       <- ssb(NSH.sams[[run]])
      ssb$run   <- run
      ssb$type  <- 'ssb'
      fbar    <- fbar(NSH.sams[[run]])
      fbar$run   <- run
      fbar$type  <- 'fbar'
      rec     <- rec(NSH.sams[[run]])
      rec$run   <- run
      rec$type  <- 'rec'
      
      df.traj <- rbind(ssb,fbar,rec)
      
      # catchabilities
      q         <- catchabilities(NSH.sams[[run]])
      q$run     <- run
      q$age[is.na(q$age)] <- "all"
      q$label   <- paste(q$fleet,q$age,sep="_")
      q$type    <- 'q'
  
      # observationv variance
      obs         <- obs.var(NSH.sams[[run]])
      obs$run     <- run
      obs$age[is.na(obs$age)] <- "all"
      obs$label   <- paste(obs$fleet,obs$age,sep="_")
      obs$type    <- 'obsVar'
  
      # uncertainties
      CV.yrs          <- ssb(NSH.sams[[run]])$year
      CV.dat          <- data.frame(year = CV.yrs,SSB=ssb(NSH.sams[[run]])$CV,
                                    Fbar=fbar(NSH.sams[[run]])$CV,Rec=rec(NSH.sams[[run]])$CV)
      CV.dat          <- CV.dat %>% pivot_longer(!year,names_to ='sub',values_to='value')
      CV.dat$run      <- run
      CV.dat$type     <- 'CV'
  
      # process variance
      mvars <- c("logSdLogFsta","logSdLogN")
      
      nSdLogFsta  <- length(which(is.element(params(NSH.sams[[run]])$name,'logSdLogFsta')))
      nSdLogN     <- length(which(is.element(params(NSH.sams[[run]])$name,'logSdLogN')))
      
      pvar        <- params(NSH.sams[[run]])[is.element(params(NSH.sams[[run]])$name,mvars),]
      pvar$run    <- run
      pvar$type   <- 'procVar'
      pvar$dummy  <- c(1:nSdLogFsta,1:nSdLogN)#1:dim(pvar)[1]
      pvar$lbnd   <- with(pvar, value - 2 * std.dev) 
      pvar$ubnd   <- with(pvar, value + 2 * std.dev)
      pvar$value  <- exp(pvar$value)
      pvar$lbnd   <- exp(pvar$lbnd)
      pvar$ubnd   <- exp(pvar$ubnd)
      pvar$paramName    <- paste(pvar$name,pvar$dummy,sep="_")
      pvar$paramType    <- gsub("log","",pvar$name)
      
      pvar <- pvar %>% select(-dummy,-name)
      
      # AIC
      df.fit      <- as.data.frame(AIC(NSH.sams[[run]]))
      df.fit$nlog <- as.data.frame(nlogl(NSH.sams[[run]]))
      df.fit$run  <- run
      colnames(df.fit) <- c('AIC','nlogl','run')
      
      if(flagFirst){
        df.all.fit <- df.fit
        df.all.pvar <- pvar
        df.all.CV   <- CV.dat
        df.all.obs  <- obs
        df.all.q    <- q
        df.all.traj <- df.traj
        
        flagFirst <- FALSE
      }else{
        df.all.fit    <- rbind(df.all.fit,df.fit)
        df.all.pvar   <- rbind(df.all.pvar,pvar)
        df.all.CV     <- rbind(df.all.CV,CV.dat)
        df.all.obs    <- rbind(df.all.obs,obs)
        df.all.q      <- rbind(df.all.q,q)
        df.all.traj   <- rbind(df.all.traj,df.traj)
      }
    }
  }
  
  return(list(df.fit = df.all.fit,
              df.q = df.all.q,
              df.obs = df.all.obs,
              df.CV = df.all.CV,
              df.pvar = df.all.pvar,
              df.traj = df.all.traj))
}

comp_sam_outputs <- function(NSH.sams){
  runNames <- names(NSH.sams)
  
  flagFirst <- TRUE
  for(run in runNames){
    if(dims(NSH.sams[[run]])$year != 1){
      # harvest
      harvest       <- as.data.frame(NSH.sams[[run]]@harvest)
      harvest$run   <- run
      harvest$type  <- 'harvest'
      # harvest
      N       <- as.data.frame(NSH.sams[[run]]@stock.n)
      N$run   <- run
      N$type  <- 'N'

      if(flagFirst){
        df.all.harvest  <- harvest
        df.all.N        <- N
        
        flagFirst <- FALSE
      }else{
        df.all.harvest  <- rbind(df.all.harvest,harvest)
        df.all.N        <- rbind(df.all.N,N)
      }
    }
  }
  
  return(list(df.harvest = df.all.harvest,
              df.N = df.all.N))
}

cor_sam_harvest <- function(NSH.sams){
  runNames <- names(NSH.sams)
  
  cor.F <- list()
  for(run in runNames){
    cor.F[[run]]    <- as.data.frame(FLCohort(NSH.sams[[run]]@harvest))
    ageUnique       <- unique(cor.F[[run]]$age)
    cor.F[[run]]    <- cor.F[[run]] %>% pivot_wider(names_from = age, values_from = data)
    cor.F[[run]]    <- log10(cor.F[[run]][,ac(ageUnique)])
  }
  
  return(cor.F)
}

span <- 9
ref.year <- 2019
type  <- 'ssb'
ref   <- '2019'

mohnRho_generic <- function(NSH.sams,ref,type){
  
  if(type == 'fbar'){
    retro  <- lapply(NSH.sams,fbar)
  }else if(type == 'ssb'){
    retro  <- lapply(NSH.sams,ssb)
  }else if(type == 'rec'){
    retro  <- lapply(NSH.sams,rec)
  }
  
  
  rho         <- data.frame(rho=NA, runs=names(retro))
  
  termFs      <- do.call(rbind,lapply(as.list(names(retro)),
                                      function(x){
                                        return(retro[[ac(x)]][which(retro[[ac(x)]]$year==max(retro[[ac(x)]]$year)),])
                                      }))
  refFs       <- retro[[ref]]
  refFs       <- refFs[which(refFs$year %in% termFs$year),]
  colnames(refFs) <- c("year","refvalue","refCV","reflbnd","refubnd")
  combFs      <- merge(refFs,termFs,by="year")
  combFs      <- combFs[order(combFs$year),]
  rhos        <- (1-combFs$value/combFs$refvalue)*100
  rho$rho     <- rhos
}

comp_sam_resi <- function(NSH.sams){
  runNames <- names(NSH.sams)
  
  flagFirst <- TRUE
  for(run in runNames){
    resi.temp       <- residuals(NSH.sams[[run]])
    resi.temp$year  <- range(NSH.sams[[1]])['minyear']+resi.temp$year-1
    range(NSH.sams[[1]])['maxyear']
    resi.temp$large <- abs(resi.temp$std.res) >= 3
    resi.temp$sign  <- ifelse(resi.temp$std.res <= 0,"N","P")
    resi.temp$run   <- run
    
    if(flagFirst){
      resi.all   <- resi.temp
      
      flagFirst <- FALSE
    }else{
      resi.all   <- rbind(resi.all,resi.temp)
    }
  }
  
  return(resi.all)
}

# ---------------------------------------------------------------------------------
# routine to fill out the output table for 95% yield intervals 
# (from Coby but organised as a subroutine)

MSY_Intervals <-function (x1.sim, interval=0.95)
{
  
  data.95 <- x1.sim$rbp
  x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
  y.95 <- data.95[data.95$variable == "Landings",]$Mean
  x.95 <- x.95[2:length(x.95)]
  y.95 <- y.95[2:length(y.95)]
  
  # Plot curve with 95% line
  ## windows(width = 10, height = 7)
  par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
       xlab = "Total catch F", ylab = "Mean landings")
  yield.p95 <- interval * max(y.95, na.rm = TRUE)
  abline(h = yield.p95, col = "blue", lty = 1)
  
  # Fit loess smoother to curve
  x.lm <- loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- predict(x.lm, newdata = lm.pred$x)
  lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")
  points(x = x1.sim$Refs["lanF","meanMSY"],
         y = predict(x.lm, newdata = x1.sim$Refs["lanF","meanMSY"]),
         pch = 16, col = "blue")
  
  # Limit fitted curve to values greater than the 95% cutoff
  lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
  fmsy.lower <- min(lm.pred.95$x)
  fmsy.upper <- max(lm.pred.95$x)
  abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
  abline(v = x1.sim$Refs["lanF","meanMSY"], lty = 1, col = "blue")
  legend(x = "bottomright", bty = "n", cex = 1.0,
         title = "F(msy)", title.col = "blue",
         legend = c(paste0("lower = ", round(fmsy.lower,3)),
                    paste0("mean = ", round(x1.sim$Refs["lanF","meanMSY"],3)),
                    paste0("upper = ", round(fmsy.upper,3))))
  
  fmsy.lower.mean <- fmsy.lower
  fmsy.upper.mean <- fmsy.upper
  landings.lower.mean <- lm.pred.95[lm.pred.95$x == fmsy.lower.mean,]$y
  landings.upper.mean <- lm.pred.95[lm.pred.95$x == fmsy.upper.mean,]$y
  
  # Repeat for 95% of yield at F(05):
  f05 <- x1.sim$Refs["catF","F05"]
  yield.f05 <- predict(x.lm, newdata = f05)
  points(f05, yield.f05, pch = 16, col = "green")
  yield.f05.95 <- interval * yield.f05
  abline(h = yield.f05.95, col = "green")
  lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
  f05.lower <- min(lm.pred.f05.95$x)
  f05.upper <- max(lm.pred.f05.95$x)
  abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
  abline(v = f05, lty = 1, col = "green")
  legend(x = "right", bty = "n", cex = 1.0,
         title = "F(5%)", title.col = "green",
         legend = c(paste0("lower = ", round(f05.lower,3)),
                    paste0("estimate = ", round(f05,3)),
                    paste0("upper = ", round(f05.upper,3))))
  
  ################################################
  # Extract yield data (landings) - median version
  
  data.95 <- x1.sim$rbp
  x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
  y.95 <- data.95[data.95$variable == "Landings",]$p50
  
  # Plot curve with 95% line
  ## windows(width = 10, height = 7)
  par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
       xlab = "Total catch F", ylab = "Median landings")
  yield.p95 <- interval * max(y.95, na.rm = TRUE)
  abline(h = yield.p95, col = "blue", lty = 1)
  
  # Fit loess smoother to curve
  x.lm <- loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- predict(x.lm, newdata = lm.pred$x)
  lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")
  
  # Find maximum of fitted curve - this will be the new median (F(msy)
  Fmsymed <- lm.pred[which.max(lm.pred$y),]$x
  Fmsymed.landings <- lm.pred[which.max(lm.pred$y),]$y
  
  # Overwrite Refs table
  x1.sim$Refs[,"medianMSY"] <- NA
  x1.sim$Refs["lanF","medianMSY"] <- Fmsymed
  x1.sim$Refs["landings","medianMSY"] <- Fmsymed.landings
  
  # Add maximum of medians to plot
  points(x = x1.sim$Refs["lanF","medianMSY"],
         y = predict(x.lm, newdata = x1.sim$Refs["lanF","medianMSY"]),
         pch = 16, col = "blue")
  
  # Limit fitted curve to values greater than the 95% cutoff
  lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
  fmsy.lower <- min(lm.pred.95$x)
  fmsy.upper <- max(lm.pred.95$x)
  abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
  abline(v = x1.sim$Refs["lanF","medianMSY"], lty = 1, col = "blue")
  legend(x = "bottomright", bty = "n", cex = 1.0,
         title = "F(msy)", title.col = "blue",
         legend = c(paste0("lower = ", round(fmsy.lower,3)),
                    paste0("median = ", round(x1.sim$Refs["lanF","medianMSY"],3)),
                    paste0("upper = ", round(fmsy.upper,3))))
  
  fmsy.lower.median <- fmsy.lower
  fmsy.upper.median <- fmsy.upper
  landings.lower.median <- lm.pred.95[lm.pred.95$x == fmsy.lower.median,]$y
  landings.upper.median <- lm.pred.95[lm.pred.95$x == fmsy.upper.median,]$y
  
  # Repeat for 95% of yield at F(05):
  f05 <- x1.sim$Refs["catF","F05"]
  yield.f05 <- predict(x.lm, newdata = f05)
  points(f05, yield.f05, pch = 16, col = "green")
  yield.f05.95 <- interval * yield.f05
  abline(h = yield.f05.95, col = "green")
  lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
  f05.lower <- min(lm.pred.f05.95$x)
  f05.upper <- max(lm.pred.f05.95$x)
  abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
  abline(v = f05, lty = 1, col = "green")
  legend(x = "right", bty = "n", cex = 1.0,
         title = "F(5%)", title.col = "green",
         legend = c(paste0("lower = ", round(f05.lower,3)),
                    paste0("estimate = ", round(f05,3)),
                    paste0("upper = ", round(f05.upper,3))))
  
  # Estimate implied SSB for each F output
  
  x.95 <- data.95[data.95$variable == "Spawning stock biomass",]$Ftarget
  b.95 <- data.95[data.95$variable == "Spawning stock biomass",]$p50
  
  # Plot curve with 95% line
  ## windows(width = 10, height = 7)
  par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  plot(x.95, b.95, ylim = c(0, max(b.95, na.rm = TRUE)),
       xlab = "Total catch F", ylab = "Median SSB")
  
  # Fit loess smoother to curve
  b.lm <- loess(b.95 ~ x.95, span = 0.2)
  b.lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                          y = rep(NA, 1000))
  b.lm.pred$y <- predict(b.lm, newdata = b.lm.pred$x)
  lines(b.lm.pred$x, b.lm.pred$y, lty = 1, col = "red")
  
  # Estimate SSB for median F(msy) and range
  b.msymed <- predict(b.lm, newdata = Fmsymed)
  b.medlower <- predict(b.lm, newdata = fmsy.lower.median)
  b.medupper <- predict(b.lm, newdata = fmsy.upper.median)
  abline(v = c(fmsy.lower.median, Fmsymed, fmsy.upper.median), col = "blue", lty = c(8,1,8))
  points(x = c(fmsy.lower.median, Fmsymed, fmsy.upper.median),
         y = c(b.medlower, b.msymed, b.medupper), col = "blue", pch = 16)
  legend(x = "topright", bty = "n", cex = 1.0,
         title = "F(msy)", title.col = "blue",
         legend = c(paste0("lower = ", round(b.medlower,0)),
                    paste0("median = ", round(b.msymed,0)),
                    paste0("upper = ", round(b.medupper,0))))
  
  # Update summary table with John's format
  
  x1.sim$Refs <- x1.sim$Refs[,!(colnames(x1.sim$Refs) %in% c("FCrash05","FCrash50"))]
  x1.sim$Refs <- cbind(x1.sim$Refs, Medlower = rep(NA,6), Meanlower = rep(NA,6),
                       Medupper = rep(NA,6), Meanupper = rep(NA,6))
  
  x1.sim$Refs["lanF","Medlower"] <- fmsy.lower.median
  x1.sim$Refs["lanF","Medupper"] <- fmsy.upper.median
  x1.sim$Refs["lanF","Meanlower"] <- fmsy.lower.mean
  x1.sim$Refs["lanF","Meanupper"] <- fmsy.upper.mean
  
  x1.sim$Refs["landings","Medlower"] <- landings.lower.median
  x1.sim$Refs["landings","Medupper"] <- landings.upper.median
  x1.sim$Refs["landings","Meanlower"] <- landings.lower.mean
  x1.sim$Refs["landings","Meanupper"] <- landings.upper.mean
  
  x1.sim$Refs["lanB","medianMSY"] <- b.msymed
  x1.sim$Refs["lanB","Medlower"] <- b.medlower
  x1.sim$Refs["lanB","Medupper"] <- b.medupper
  
  # Reference point estimates
  cat("\nReference point estimates:\n")
  return (round(x1.sim$Refs,3))
}


# -End of function-------------------------------------------------------------------------

# modified from msy::eqsr_fit to allow shifting in S-R time series
eqsr_fit_shift <- 
  function (stk, nsamp = 5000, models = c("ricker", "segreg", "bevholt"), 
            method = "Buckland", id.sr = NULL, remove.years = NULL, delta = 1.3, 
            nburn = 10000, rshift = 0) 
  {
    dms <- FLCore::dims(stk)
    rage <- dms$min
    if (rage == 0) {
      x = FLCore::stock.n(stk)[1, drop = TRUE]
    }
    else {
      x = c(FLCore::stock.n(stk)[1, -seq(rage), drop = TRUE], 
            rep(NA, rage))
    }
    if (rshift > 0){
      x = c(FLCore::stock.n(stk)[1, -seq(rshift), drop = TRUE], 
            rep(NA, rshift))
      
    } else { NULL }        
    rby <- data.frame(year = with(dms, minyear:maxyear), rec = x, 
                      ssb = FLCore::ssb(stk)[drop = TRUE], fbar = FLCore::fbar(stk)[drop = TRUE], 
                      landings = FLCore::landings(stk)[drop = TRUE], catch = FLCore::catch(stk)[drop = TRUE])
    # print(rby)
    row.names(rby) <- NULL
    rby <- rby[!is.na(rby$rec), ]
    data <- rby[, 1:3]
    if (!is.null(remove.years)) {
      data$ssb[data$year %in% remove.years] <- NA
    }
    data <- data[complete.cases(data), ]
    if (is.null(id.sr)) 
      id.sr <- FLCore::name(stk)
    method <- match.arg(method, c("Buckland", "Simmonds", "King", "Cadigan"))
    if (!is.character(models)) 
      stop("models arg should be character vector giving names of stock recruit models")
    if (method == "Buckland") {
      return(c(eqsr_Buckland(data, nsamp, models), list(stk = stk, 
                                                             rby = rby, id.sr = id.sr)))
    }
    else {
      cat("The", method, "is not ready yet!  Working on it!\n")
    }
  }

# -End of function-------------------------------------------------------------------------

#' stock recruitment function
#'
#'
#' @param data data.frame containing stock recruitment data
#' @param nsamp Number of samples
#' @param models A character vector
#' @param ... Additional arguements
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@ices.dk}
#' @export
eqsr_Buckland <- function(data, nsamp = 5000, models = c("Ricker","Segreg","Bevholt"), ...)
{
  
  ## dummy
  model <- 0
  #--------------------------------------------------------
  # Fit models
  #--------------------------------------------------------
  
  nllik <- function(param, ...) -1 * llik(param, ...)
  ndat <- nrow(data)
  fit <- lapply(1:nsamp, function(i)
  {
    sdat <- data[sample(1:ndat, replace = TRUE),]
    
    fits <- lapply(models, function(mod) stats::nlminb(initial(mod, sdat), nllik, data = sdat, model = mod, logpar = TRUE))
    
    best <- which.min(sapply(fits, "[[", "objective"))
    
    with(fits[[best]], c(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = best))
  })
  
  fit <- as.data.frame(do.call(rbind, fit))
  fit $ model <- models[fit $ model]
  
  #--------------------------------------------------------
  # get posterior distribution of estimated recruitment
  #--------------------------------------------------------
  pred <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], sort(data $ ssb))) ))
  #pred <- t(sapply(seq(nsamp), function(j) exp(get(fit $ model[j], , pos = "package:msy", mode = "function") (fit[j,], sort(data $ ssb))) ))
  
  
  #--------------------------------------------------------
  # get best fit for each model
  #--------------------------------------------------------
  fits <-
    do.call(rbind,
            lapply(models,
                   function(mod)
                     with(stats::nlminb(initial(mod, data), nllik, data = data, model = mod, logpar = TRUE),
                          data.frame(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = mod))))
  
  tmp <- plyr::ddply(fit,c("model"), plyr::summarise, n=length(model))
  tmp$prop <- tmp$n/sum(tmp$n)
  fits <- plyr::join(fits,tmp,by="model")
  
  dimnames(pred) <- list(model=fit$model,ssb=data$ssb)
  
  return(list(sr.sto = fit, sr.det = fits, pRec = pred))
}

# -End of function-------------------------------------------------------------------------
