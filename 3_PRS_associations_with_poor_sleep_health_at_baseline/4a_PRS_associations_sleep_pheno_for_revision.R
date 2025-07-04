library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)


dat <- read.csv("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/20250514_cleaned_DM_PRS_data_w_covars.csv")

#convert from integer to a factor variable
dat$EDUCATION_C3 <- factor(dat$EDUCATION_C3)
dat$EDUCATION_C3 <- factor(dat$MED_STATIN)
dat$EDUCATION_C3 <- factor(dat$INCOME_C3)


outcomes <- c("Mild_OSA", "Mod_severe_OSA", "Mild_severe_OSA")


PRSs <- c("PRSstd_sum_std",
          "PRSstd_gap_std",
          "PRSstd_mgb_std")

#First, replace BMI with WHR as per reviewers comment
#"The authors used BMI as an obesity metric. 
# Considering the role of visceral fat in mediating the results. 
# Would replacing BMI with a measure of central obesity such as abdominal circumference 
# or waist-to-hip ratio potentially provide more robust evidence on the strength of association in the model?"

res <- c()

for (i in 1:length(outcomes)){
  outcome <- outcomes[i]
  dat[[outcome]] <- ifelse(dat[[outcome]] == "Yes", 1, 0)
  dat$keep <- ifelse( !is.na(dat$PRS_sum),1, 0)
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    cur_prs <- PRSs[j]
    
    
    mod <- svyglm(as.formula(paste(outcome, "~ Age + Gender + WAIST_HIP + Center +", 
                                   "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +", 
                                   cur_prs)), 
                  design = survey_obj, 
                  family = "quasibinomial")
    
    res <- rbind(res, data.frame(
      outcome = outcome,
      n = nrow(model.matrix(mod)),
      N_case = table(mod$y)[2],
      N_control = table(mod$y)[1],
      PRS = cur_prs,
      log_OR = summary(mod)$coef[cur_prs,"Estimate"],
      log_PR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
      log_OR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"]))
    
    
    
    
  }
}

write.csv(res, file = "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Code/Code_for_GitHub/Results/20250514_assoc_PRS_sleep_pheno_w_WHR.csv", row.names = FALSE)


#these are the variables to be adjusted for one at a time as requested by one of the reviewers
vars_adj_for <- c("EDUCATION_C3", "MED_STATIN", "INCOME_C3",
                  "GPAQ_TOTAL_VIG", "GPAQ_TOTAL_MET", "BKGRD1_C7")

res2 <- c()

# Use baseline data to perform association of diabetes PRS with sleep disturbances
for (i in 1:length(outcomes)){
  
  outcome <- outcomes[i]
  
  dat[[outcome]] <- ifelse(dat[[outcome]] == "Yes", 1, 0)
  dat$keep <- ifelse( !is.na(dat$PRS_sum),1, 0)
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    
    cur_prs <- PRSs[j]
    
    for (k in 1:length(vars_adj_for)){
      
      cur_var <- vars_adj_for[k]
      
      mod <- svyglm(as.formula(paste(outcome, "~ Age + Gender + BMI + Center +", 
                                     "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +", cur_var, "+", cur_prs)), 
                    design = survey_obj, 
                    family = "quasibinomial")
      
  
      
      res2 <- rbind(res2, data.frame(
        outcome = outcome,
        var_adj = cur_var,
        n = nrow(model.matrix(mod)),
        N_case = table(mod$y)[2],
        N_control = table(mod$y)[1],
        PRS = cur_prs,
        log_OR = summary(mod)$coef[cur_prs,"Estimate"],
        log_PR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
        log_OR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"]))
      
      
      
    }
    
  }
}





#write.csv(res2, file = "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Code/Code_for_GitHub/Results/20250514_assoc_PRS_sleep_pheno_w_additional_confounders.csv", row.names = FALSE)





