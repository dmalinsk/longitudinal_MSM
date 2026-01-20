## uses the sl3 R package (ver 1.4.5)
## install with: remotes::install_github("tlverse/sl3")
library(sl3)

sl_density <- function(data, exposure, covariates){

  ok <- complete.cases(data[,c(exposure,covariates)])
  data <- data[ok,]
  
  task <- make_sl3_Task(
    data = data,
    outcome = exposure,
    covariates = covariates
  )

  hese_glm_glm_lrnr <- Lrnr_density_semiparametric$new(name="GLM_GLM",
                                                     mean_learner = Lrnr_glm$new(),
                                                     var_learner = Lrnr_glm$new()
  )

  hese_rf_rf_lrnr <- Lrnr_density_semiparametric$new(name="RF_RF",
                                                   mean_learner = Lrnr_ranger$new(),
                                                   var_learner = Lrnr_ranger$new()
  )

  # hese_hal_lrnr <- Lrnr_density_semiparametric$new(name="HAL",
  #   mean_learner = Lrnr_hal9001$new(),
  #   var_learner = Lrnr_hal9001$new()
  # )

  hese_xgboost_lrnr <- Lrnr_density_semiparametric$new(name="xgboost",
                                                     mean_learner = Lrnr_xgboost$new(),
                                                     var_learner = Lrnr_xgboost$new()
  )

  hese_gam_lrnr <- Lrnr_density_semiparametric$new(name="GAM",
                                                 mean_learner = Lrnr_gam$new(),
                                                 var_learner = Lrnr_gam$new()
  )

  hese_lasso_lrnr <- Lrnr_density_semiparametric$new(name="lasso",
                                                   mean_learner = Lrnr_glmnet$new(alpha=1),
                                                   var_learner = Lrnr_glmnet$new(alpha=1)
  )

  hese_svm_lrnr <- Lrnr_density_semiparametric$new(name="SVM",
                                                 mean_learner = Lrnr_svm$new(),
                                                 var_learner = Lrnr_svm$new()
  )

  # SL for the conditional treatment density
  sl_dens_lrnr <- Lrnr_sl$new(
    learners = list(hese_glm_glm_lrnr, hese_rf_rf_lrnr,
                  hese_xgboost_lrnr, hese_gam_lrnr, hese_lasso_lrnr, hese_svm_lrnr),
    metalearner = Lrnr_solnp_density$new()
  )

  sl_fit <- sl_dens_lrnr$train(task = task)
  pred <- sl_fit$predict(task = task)
  
  w <- rep(NA,length(data[,exposure]))
  w[ok] <- pred$likelihood

  return(w)
  
}