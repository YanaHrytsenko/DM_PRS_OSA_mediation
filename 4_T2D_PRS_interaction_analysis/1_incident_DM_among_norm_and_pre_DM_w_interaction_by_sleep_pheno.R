library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)


dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")




PRSs <- c("PRSstd_sum_std",
          "PRSstd_gap_std",
          "PRSstd_mgb_std")

sleep_exposures <-  c("Mild_OSA", "Mod_severe_OSA", "Short_sleep", "Long_sleep", "EDS", "Insomnia")


# for each strata of interest, compute the association of DM PRS with incident diabetes. 
# here we would limit the analysis to individuals who participated in visit 2, and use visit 2 weights. 
dat <- dat[which(!is.na(dat$WEIGHT_NORM_OVERALL_V2)),]


# create a column to indicate inds w DM at visit 2
dat$inc_diab <- 0
dat$inc_diab[which(dat$Visit2_diabetes == "Diabetic")] <- 1
dat$inc_diab[which(is.na(dat$Visit2_diabetes == "Diabetic"))] <- NA 



#Primary incident analysis: normal and pre-diabetic at baseline
res <- c()

dat$keep <- ifelse( !is.na(dat$PRS_sum) &
                      (dat$Baseline_diabetes == "Normoglycemic" | dat$Baseline_diabetes == "Hyperglycemic"),1, 0)

survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_V2, data=dat)
survey_obj <- subset(survey_obj, keep == 1)



for (j in 1:length(PRSs)){
  #for each PRS of interest
  cur_prs <- PRSs[j]
  
  #for each sleep phenotype
  for (k in 1:length(sleep_exposures)){
    
    mod <- svyglm(as.formula(paste("inc_diab ~ Age + Gender + BMI + Center +",
                                   "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + ",
                                   paste(paste(sleep_exposures[k], cur_prs, sep = "*"), collapse = "+"), "+",
                                   cur_prs)),
                  
                  offset = log(Time_between_visits), 
                  design = survey_obj, 
                  family = quasipoisson(link = "log"))
    
    
    exposures <- setdiff(names(mod$coef), c("(Intercept)", "Age", 'GenderMale', 
                                            "BMI", "CenterChicago", "CenterMiami",
                                            "CenterSan Diego", "PC_1", 
                                            "PC_2", "PC_3", "PC_4", "PC_5"))
    
    res <- rbind(res, data.frame(
      model = "sleep_interactions",
      n = nrow(model.matrix(mod)),
      PRS = cur_prs,
      Exposure = exposures,
      log_IRR = summary(mod)$coef[exposures,"Estimate"],
      log_IRR_SE = summary(mod)$coef[exposures,"Std. Error"],
      log_IRR_p = summary(mod)$coef[exposures,"Pr(>|t|)"])
    )
  }
  
}



write.csv(res, file = "/Results/assoc_PRS_incident_DM_norm_and_pre_DM_sleep_adjusted_w_interaction.csv")
