library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)
library(data.table)

dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")

prs_BMI_adj <- fread("/Data/EUR_BMI_adjusted_diabetes_GWAS/20240404_HCHS_BMIadj_T2D_PRS_CS.txt", header = T)

prs_BMI_adj <- prs_BMI_adj[, c("IID", "std_score_BMIadj_T2D_PRS_CS")]

colnames(prs_BMI_adj) <- c("SUBJECT_ID", "BMIadj_T2D_PRS")

dat <- left_join(dat, prs_BMI_adj, by = "SUBJECT_ID")



outcomes <- c("Mild_severe_OSA")

PRSs <- c("PRSstd_sum_std",
          "PRSstd_gap_std",
          "PRSstd_mgb_std", 
          "BMIadj_T2D_PRS")


# Use baseline data to perform association of diabetes PRS with sleep disturbances
res <- c()

for (i in 1:length(outcomes)){
  
  outcome <- outcomes[i]
  dat[[outcome]] <- ifelse(dat[[outcome]] == "Yes", 1, 0)
  dat$keep <- ifelse( !is.na(dat$PRS_sum),1, 0)
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    
    cur_prs <- PRSs[j]
    
    
    mod <- svyglm(as.formula(paste(outcome, "~ Age + Gender + BMI + Center +", 
                                   "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +", 
                                   cur_prs)), 
                  design = survey_obj, 
                  family = "quasibinomial")
    
    res <- rbind(res, data.frame(
      outcome = outcome,
      n= nrow(model.matrix(mod)),
      N_case = table(mod$y)[2],
      N_control = table(mod$y)[1],
      PRS = cur_prs,
      log_OR = summary(mod)$coef[cur_prs,"Estimate"],
      log_PR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
      log_OR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"]))
    
    
    
  }
}



write.csv(res, file = "/Results/assoc_PRS_OSA.csv")





