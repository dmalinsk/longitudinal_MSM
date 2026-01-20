library(haven)
library(lme4)

emph <- data.frame() ## load data here

attach(emph)
sample = wrong_CT==0 & aquilion_CT==0 ## this is because the JAMA analysis excludes data with wrong_CT or aquilon_CT=1
detach(emph)
emph <- emph[sample,]
rm(sample)

emph$site1c <- as.factor(emph$site1c)
emph$site <- as.factor(emph$site)
emph$edu <- as.factor(emph$edu)
emph$exercm_cat <- as.factor(emph$exercm_cat)

### Some NOx and PM2.5 values are NA when they shouldn't be. Carrying-over first value. ###
for(i in 2:dim(emph)[1]){
  if( emph[i,"idno"]==emph[i-1,"idno"] ) {
    # emph[i,"o3_base"] <- emph[i-1,"o3_base"] ## this is fixed after data update
    emph[i,"nox_base"] <- emph[i-1,"nox_base"]
    emph[i,"pm25_base"] <- emph[i-1,"pm25_base"] ## for NOx and PM2.5, sometimes row i is NA when i-1 is not-NA
  }
}
###


## begin replication of JAMA result
#######################

## multiply coefs by 100 for NOx, 20 for PM2.5, and 30 for ozone

# mod.jama <- lmer(B_emph950_out ~ timefactor + o3_base + timefactor:o3_follow
#              # time-varying vars
#              + New_smkst:timefactor + Cigsday16_t:timefactor + shndsmk1:timefactor
#              + Htcm16:timefactor + Wtlb16:timefactor + Bmi16:timefactor
#              #+ temp_follow:timefactor ## JAMA analysis had temp
#              + nsesindex:site
#              # baseline vars
#              + site1c + Age1c + Black + Chinese + Hispanic + Gender1 + htcm16_c + wtlb16_c + bmi16_c #+ temp_base ## JAMA analysis had temp, removing now
#              + site:timefactor + Age1c:timefactor + Black:timefactor + Chinese:timefactor + Hispanic:timefactor + Gender1:timefactor
#              + edu*timefactor + income*timefactor + exercm_cat*timefactor + employ*timefactor + nsesindex
#              # scanner vars
#              + Scanner_EBT + Scanner_JHUVZ4 + Scanner_eMDCTnoJH + Scanner_Force
#              + Scanner_GE64 + Scanner_MNLA64 + Vxsz_11 + Vxsz_14
#              # misc (smoking and scanner)
#              + PY1_A + PY1_A:timefactor + Miss_PY1_A + Miss_PY1_A:timefactor
#              + new_smkst_c + cigsday16_t_c + shndsmk1_c + Wt220 + High_BMI5 + Low_BMI5
#              # random intercept and slope
#              + (1 + timefactor | idno), data = emph, REML = FALSE, control = lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e6))) # , corr=corCAR1(form=~timefactor|idno))
# summary(mod.jama)
# 
# confint(mod.jama,"timefactor:o3_follow")

## end replication
#######################


## making a new exposure variable for ozone that combines baseline and follow up values
## also a variable which captures cumulative sum of exposure
o3_exp <- c()
nox_exp <- c()
pm25_exp <- c()
o3_sum <- c()
nox_sum <- c()
pm25_sum <- c()

for(i in 1:dim(emph)[1]){
  o3_exp[i] <- ifelse(emph[i,"o3_follow"]==0,emph[i,"o3_base"],emph[i,"o3_follow"])
  nox_exp[i] <- ifelse(emph[i,"nox_follow"]==0,emph[i,"nox_base"],emph[i,"nox_follow"])
  pm25_exp[i] <- ifelse(emph[i,"pm25_follow"]==0,emph[i,"pm25_base"],emph[i,"pm25_follow"])
  
  o3_sum[i] <- ifelse(emph[i,"o3_follow"]==0,emph[i,"o3_base"],emph[i,"o3_base"]+emph[i,"o3_follow"])
  nox_sum[i] <- ifelse(emph[i,"nox_follow"]==0,emph[i,"nox_base"],emph[i,"nox_base"]+emph[i,"nox_follow"])
  pm25_sum[i] <- ifelse(emph[i,"pm25_follow"]==0,emph[i,"pm25_base"],emph[i,"pm25_base"]+emph[i,"pm25_follow"])
}

emph <- as.data.frame(emph)

emph$o3_exp <- as.numeric(o3_exp)
emph$nox_exp <- as.numeric(nox_exp)
emph$pm25_exp <- as.numeric(pm25_exp)
emph$combo_exp <- as.numeric(emph$o3_exp + emph$nox_exp + emph$pm25_exp)

emph$o3_sum <- as.numeric(o3_sum)
emph$nox_sum <- as.numeric(nox_sum)
emph$pm25_sum <- as.numeric(pm25_sum)
emph$combo_sum <- as.numeric(emph$o3_sum + emph$nox_sum + emph$pm25_sum)

rm(o3_exp)
rm(nox_exp)
rm(pm25_exp)
rm(o3_sum)
rm(nox_sum)
rm(pm25_sum)
## end variable construction


## calculate weights
#######################
statcovariates <- c("site1c", "Age1c", "Black", "Chinese", "Hispanic", "Gender1", 
                    "htcm16_c", "wtlb16_c", "bmi16_c", 
                    "PY1_A",  "Miss_PY1_A",  "new_smkst_c", "cigsday16_t_c", "shndsmk1_c", 
                    "edu", "income", "exercm_cat", "employ", "nsesindex"
                    , "pm25_base","nox_base"
                    )
tvcovariates <- c( 
                  "Bmi16",
                  "Cigsday16_t", "New_smkst", "shndsmk1"
                  , "pm25_follow", "nox_follow"
                  )

source("ipw.R")

weights <-  ipw(data=emph, 
                exposure="o3_exp", 
                timevar="ct_exam", 
                id="idno", 
                quantiles=5, 
                numvisits=6,
                statcovariates=statcovariates,
                tvcovariates=tvcovariates,
                outcome="B_emph950_out",
                timefactor="timefactor",
                trunc=TRUE,
                trunc_level=0.01,
                statinteractions = c("nsesindex:site1c"),
                mode = "binning") # two modes: "binning" or "sl_density"

summary(weights$sw)
emph$sw <- weights$sw

mod <- lmer(B_emph950_out ~ timefactor + o3_base + timefactor:o3_follow
                 # scanner vars
                 + Scanner_EBT + Scanner_JHUVZ4 + Scanner_eMDCTnoJH + Scanner_Force
                 + Scanner_GE64 + Scanner_MNLA64 + Vxsz_11 + Vxsz_14
                 + Wt220 + High_BMI5 + Low_BMI5
                 # random intercept and slope
                 + (1 + timefactor | idno), 
                 data = emph, REML = TRUE, weights = emph$sw, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

summary(mod)

###########################################################################