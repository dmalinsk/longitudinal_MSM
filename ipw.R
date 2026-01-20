library(dplyr)
library(nnet)
source("sl_density.R")

ipw <- function(
  exposure,
  id = NULL,
  timevar = NULL,
  numvisits = NULL,
  data,
  statcovariates = NULL,
  tvcovariates = NULL,
  outcome = NULL,
  trunc = NULL,
  quantiles = 10,
  timefactor = NULL,
  trunc_level = NULL,
  statinteractions = NULL,
  mode = "binning"
)

  {
  tempcall <- match.call()
  
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("id" %in% names(tempcall))) stop("No patient id specified")
  if (!("timevar" %in% names(tempcall))) stop("No timevar specified")
  if (!("numvisits" %in% names(tempcall))) stop("No numvisits specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("statcovariates" %in% names(tempcall))) stop("No static covariates specified")
  if (!("tvcovariates" %in% names(tempcall))) stop("No time-varying covariates specified")
  if (!("outcome" %in% names(tempcall))) stop("No outcome specified")
  if (!("timefactor" %in% names(tempcall))) stop("No timefactor specified")
  
  #if (!is.null(tempcall$trunc)) {if(tempcall$trunc < 0 | tempcall$trunc > 0.5) stop("Invalid truncation percentage specified (0-0.5)")}
  
  tempdat <- data
  # make wide format
  tempdat <- as.data.frame(tempdat[,c(id, timevar, timefactor, statcovariates, tvcovariates, exposure, outcome)])
  dataW <- reshape(tempdat, v.names = c(timefactor, tvcovariates, exposure, outcome), timevar = timevar, idvar = id, direction = "wide")
    
  
  adj <- function(x){
    x[which(is.na(x))] <- 1
    return(x)
  }

  ## bin the exposure into quantiles. Default is 10.

  for (i in 1:numvisits) {
    dataW[paste0("quant.",i)] = as.numeric(cut(dataW[,paste(exposure, i, sep = ".")], quantile(dataW[,paste(exposure, i, sep = ".")], probs=0:quantiles/quantiles, na.rm = T), include.lowest = T))
    }
  
  # calculate weights
  
  # create model for first visit

  if(mode == "binning"){
  
    statcovariates <- c(statcovariates,statinteractions) ## to force specified interaction
    
    mod.1 <- multinom(paste("quant.1 ~ ",paste(statcovariates, collapse="+"),sep = ""), data = dataW)
  
    propensA <- predict(mod.1,type="probs",newdata = dataW) ## gives a N x 10 matrix, have to select observed propens
    pA1 <- NULL
    for(i in 1:nrow(dataW)){
      pA1 <- rbind(pA1,propensA[i,dataW$quant.1[i]])
    }
  
    w.1 <- as.numeric((1/quantiles)/pA1)
  } #end if(mode == "binning")
  
  if(mode == "sl_density"){
    den.1 <- sl_density(data=dataW,exposure=paste(exposure, 1, sep = "."),covariates=statcovariates)
    exposure.1 <- paste(exposure, 1, sep = ".")
    complete_exp <- dataW[complete.cases(dataW[,c(exposure.1,statcovariates)]), exposure.1]
    num.1 <- rep(NA,length(dataW[,exposure.1]))
    num.1[complete.cases(dataW[,c(exposure.1,statcovariates)])] <- approx(density(complete_exp)$x, density(complete_exp)$y, xout = complete_exp)$y
    w.1 <- num.1/den.1
  } #end if(mode == "sl_density")
  
  ## cat("weights dist at baseline =", summary(w.1),"\n")
  
  ### truncating weights before taking the product
  if (trunc == TRUE){
    w.1.trunc <- w.1
    w.1.trunc[w.1 <= quantile(w.1, 0+trunc_level,na.rm=T)] <- quantile(w.1, 0+trunc_level,na.rm=T)
    w.1.trunc[w.1 >  quantile(w.1, 1-trunc_level,na.rm=T)] <- quantile(w.1, 1-trunc_level,na.rm=T)
    w.1 <- w.1.trunc / mean(w.1.trunc,na.rm=T) ### makes mean of weights =1
    ## cat("truncated weights dist at baseline = ", summary(w.1.trunc),"\n")
  }
  ###
  
  ipw.1 <- w.1
  prevwts <- adj(w.1)
  dataW$ipw.1 <- ipw.1
  
  models <- list()
  weights <- list()
  baseline_exp <- paste(exposure,".1", sep="")
  baseline_outcome <- paste(outcome,".1", sep="")
  num.models <- list()
  
  ## create models for all visits after 1 and calculate weights
  
  for (visit in 2:numvisits) {
    
    if(mode == "binning"){
    
      visit_quant <- paste("quant.", visit, sep = "")
      visit_tvcovars <- paste(tvcovariates,".", visit, "*", timefactor,".", visit, sep="")
      
      models[[visit - 1]] <- multinom(paste("quant.", visit,"~", 
                                          paste(statcovariates, collapse="+"), 
                                          paste("+", visit_tvcovars, collapse="+"), 
                                          #paste("+", baseline_exp, sep = ""),
                                          #paste("+", baseline_exp, "*", timefactor, ".", visit, sep = ""),
                                          paste("+", baseline_outcome, sep = ""),
                                          sep=""),
                                    data = dataW)
    
      propensA <- predict(models[[visit - 1]],type="probs",newdata = dataW) ## gives a N x q matrix, have to select observed propens
      quant_dat <- dataW[visit_quant] ## create intermediate variable to access quantile values in loop
      pA <- NULL
      for (i in 1:nrow(dataW)) {
        pA <- rbind(pA,propensA[i,quant_dat[i,]])
      }
    
      num.models[[visit - 1]] <- multinom(paste("quant.", visit,"~", 
                                          paste(statcovariates, collapse="+"), 
                                          sep=""),
                                    data = dataW)
    
      propensA.num <- predict(num.models[[visit - 1]],type="probs",newdata = dataW) ## gives a N x q matrix, have to select observed propens
      pA.num <- NULL
      for (i in 1:nrow(dataW)) {
        pA.num <- rbind(pA.num,propensA.num[i,quant_dat[i,]])
      }
    
      weights[[visit - 1]] <- as.numeric(pA.num/pA)
    
    } #end if(mode == "binning")
    
    if(mode == "sl_density"){
      
      visit_tvcovariates <- paste(tvcovariates,".",visit,sep="")
      
      den.v <- sl_density(data = dataW, exposure=paste(exposure, visit, sep = "."),
                        covariates=c(statcovariates,visit_tvcovariates,baseline_outcome,paste(timefactor,".",visit,sep="")))
      num.v <- sl_density(data = dataW, exposure=paste(exposure, visit, sep = "."), 
                        covariates=c(statcovariates))
      weights[[visit - 1]] <- num.v/den.v
    } #end if(mod == "sl_density)
    
    wt <- weights[[visit-1]]
    
    ## cat("weights dist for visit ", visit, " = ", summary(wt),"\n")
    
    ### truncating weights before taking the product
    if (trunc == TRUE){
      wt.trunc <- wt
      wt.trunc[wt <= quantile(wt, 0+trunc_level,na.rm=T)] <- quantile(wt, 0+trunc_level,na.rm=T)
      wt.trunc[wt >  quantile(wt, 1-trunc_level,na.rm=T)] <- quantile(wt, 1-trunc_level,na.rm=T)
      wt <- wt.trunc / mean(wt.trunc, na.rm=T) ### makes mean of weights =1
      ## cat("truncated weights dist for visit ", visit, " = ", summary(wt.trunc),"\n")
    }
    ###
    assign(paste("ipw",visit, sep="."), prevwts*wt)
    dataW[paste("ipw",visit, sep=".")] <- eval(as.name(paste("ipw",visit, sep=".")))
    prevwts <- prevwts*adj(wt)
  }
  
  sw <- rep(NA,length(tempdat[,id]))
  
  for (item in 1:numvisits) {
    for (i in 1:length(dataW[,id])){
      if (!is.na(dataW[i,paste("ipw",item,sep=".")])){
        sw[which(tempdat[,id] == dataW[,id][i] & tempdat[,timevar] == item)] <- dataW[i,paste("ipw",item,sep=".")]
      }
    }
  }
  
  
  weights <- data.frame(sw)
  
  
  return(weights)
}