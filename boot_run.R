## THIS BIT OF CODE BOOTSTRAPS THE MAIN ANALYSIS FOR EACH POLLUTANT (multipollutant model) ##
## ASSUMES YOU HAVE THE PRE-PROCESSED DATA "emph" (from analysis.R file) AS INPUT ##
## MAY BE RUN WITH METHOD = "binning" OR "sl_density" ##

source("ipw.R")

boot_sample <- function(x) {
  temp <- emph[emph$idno == x,]
  assign("newid", value = newid + 1, envir = .GlobalEnv)
  # newid <- newid + 1
  temp$newid <- newid
  return(temp)
}

niter <- 1000 ## (Note: SuperLearner takes a long time, so niter=200 for "sl_density" vs. 1000 for "binning") 

## ozone ##

set.seed(3)

coef_list <- c()

for (iter in 1:niter) { ########TEMP#########
  newid <- 0
  samp <- sample(unique(emph$idno), size=7071, replace=TRUE)
  
  resample <- do.call(rbind,lapply(samp, boot_sample))
  
  cat("Estimating weights for ozone, iteration # = ", iter, "\n")
  
  weights <- ipw(data=resample, 
                    exposure="o3_exp", 
                    timevar="ct_exam", 
                    id="newid", 
                    quantiles=5, 
                    numvisits=6,
                    statcovariates <- c("site1c", "Age1c", "Black", "Chinese", "Hispanic", "Gender1", 
                                        "htcm16_c", "wtlb16_c", "bmi16_c", 
                                        "PY1_A",  "Miss_PY1_A",  "new_smkst_c", "cigsday16_t_c", "shndsmk1_c", 
                                        "edu", "income", "exercm_cat", "employ", "nsesindex"
                                        , "pm25_base", "nox_base"
                    ),
                    tvcovariates <- c( 
                      "Bmi16", 
                      "Cigsday16_t", "New_smkst", "shndsmk1"
                      , "pm25_follow", "nox_follow"
                    ),
                    outcome="B_emph950_out",
                    timefactor="timefactor",
                    trunc=TRUE,
                    trunc_level=0.01,
                    statinteractions = c("nsesindex:site1c"),
                    mode = "binning")
  
  
  resample$sw <- weights$sw
  
  mod <- lmer(B_emph950_out ~ timefactor + o3_base + timefactor:o3_follow
              # scanner vars
              + Scanner_EBT + Scanner_JHUVZ4 + Scanner_eMDCTnoJH + Scanner_Force
              + Scanner_GE64 + Scanner_MNLA64 + Vxsz_11 + Vxsz_14
              + Wt220 + High_BMI5 + Low_BMI5
              # random intercept and slope
              + (1 + timefactor | newid), 
              data = resample, REML = TRUE, weights = resample$sw, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  coef_list[iter] <- coef(summary(mod))[15] ## don't forget to check this! relevant coef changes w/ model specification
}

res.o3 <- quantile(coef_list, c(0.025, 0.5, 0.975))
sd.o3 <- sqrt(var(coef_list))
coef_list.o3 <- coef_list

## end ozone ##

## nox ##

set.seed(3)

coef_list <- c()

for (iter in 1:niter) {
  newid <- 0
  samp <- sample(unique(emph$idno), size=7071, replace=TRUE)
  
  resample <- do.call(rbind,lapply(samp, boot_sample))
  
  cat("Estimating weights for nox, iteration # = ", iter, "\n")
  
  weights <- ipw(data=resample, 
                    exposure="nox_exp", 
                    timevar="ct_exam", 
                    id="newid", 
                    quantiles=5, 
                    numvisits=6,
                    statcovariates <- c("site1c", "Age1c", "Black", "Chinese", "Hispanic", "Gender1", 
                                        "htcm16_c", "wtlb16_c", "bmi16_c", 
                                        "PY1_A",  "Miss_PY1_A",  "new_smkst_c", "cigsday16_t_c", "shndsmk1_c", 
                                        "edu", "income", "exercm_cat", "employ", "nsesindex"
                                        , "pm25_base", "o3_base"
                    ),
                    tvcovariates <- c( 
                      "Bmi16", 
                      "Cigsday16_t", "New_smkst", "shndsmk1"
                      , "pm25_follow", "o3_follow"
                    ),
                    outcome="B_emph950_out",
                    timefactor="timefactor",
                    trunc=TRUE,
                    trunc_level=0.01,
                    statinteractions = c("nsesindex:site1c"),
                    mode = "binning")
  
  
  resample$sw <- weights$sw
  
  mod <- lmer(B_emph950_out ~ timefactor + nox_base + timefactor:nox_follow
              # scanner vars
              + Scanner_EBT + Scanner_JHUVZ4 + Scanner_eMDCTnoJH + Scanner_Force
              + Scanner_GE64 + Scanner_MNLA64 + Vxsz_11 + Vxsz_14
              + Wt220 + High_BMI5 + Low_BMI5
              # random intercept and slope
              + (1 + timefactor | newid), 
              data = resample, REML = TRUE, weights = resample$sw, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  coef_list[iter] <- coef(summary(mod))[15] ## don't forget to check this! relevant coef changes w/ model specification
}

res.nox <- quantile(coef_list, c(0.025, 0.5, 0.975))
sd.nox <- sqrt(var(coef_list))
coef_list.nox <- coef_list

## end nox ##

## pm25 ##
set.seed(3)

coef_list <- c()

for (iter in 1:niter) {
  newid <- 0
  samp <- sample(unique(emph$idno), size=7071, replace=TRUE) ## full sample size 7071, 6364 is 90%

  resample <- do.call(rbind,lapply(samp, boot_sample))

  cat("Estimating weights for pm25, iteration # = ", iter, "\n")

  weights <- ipw(data=resample,
                    exposure="pm25_exp",
                    timevar="ct_exam",
                    id="newid",
                    quantiles=5,
                    numvisits=6,
                    statcovariates <- c("site1c", "Age1c", "Black", "Chinese", "Hispanic", "Gender1",
                                        "htcm16_c", "wtlb16_c", "bmi16_c",
                                        "PY1_A",  "Miss_PY1_A",  "new_smkst_c", "cigsday16_t_c", "shndsmk1_c",
                                        "edu", "income", "exercm_cat", "employ", "nsesindex"
                                        , "o3_base", "nox_base"
                    ),
                    tvcovariates <- c(#"Htcm16", "Wtlb16",
                      "Bmi16",
                      "Cigsday16_t", "New_smkst", "shndsmk1"
                      , "o3_follow", "nox_follow"
                    ),
                    outcome="B_emph950_out",
                    timefactor="timefactor",
                    trunc=TRUE,
                    trunc_level=0.01,
                    statinteractions = c("nsesindex:site1c"),
                    mode = "binning")


  resample$sw <- weights$sw

  mod <- lmer(B_emph950_out ~ timefactor + pm25_base + timefactor:pm25_follow
              # scanner vars
              + Scanner_EBT + Scanner_JHUVZ4 + Scanner_eMDCTnoJH + Scanner_Force
              + Scanner_GE64 + Scanner_MNLA64 + Vxsz_11 + Vxsz_14
              + Wt220 + High_BMI5 + Low_BMI5
              # random intercept and slope
              + (1 + timefactor | newid),
              data = resample, REML = TRUE,
              weights = resample$sw,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

  coef_list[iter] <- coef(summary(mod))[15] ## don't forget to check this! relevant coef changes w/ model specification
}

res.pm25 <- quantile(coef_list, c(0.025, 0.5, 0.975))
sd.pm25 <- sqrt(var(coef_list))
coef_list.pm25 <- coef_list

## end pm25 ##