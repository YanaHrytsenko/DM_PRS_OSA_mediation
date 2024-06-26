library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)
library(data.table)

dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")

OSA_PRS_bmi_UNadj <- fread("/Data/OSA_PRS/bmi_UNadj/OSA_bmi_UNadj.txt", header = T)

OSA_PRS_bmi_UNadj <- OSA_PRS_bmi_UNadj[,-c("FID")]

colnames(OSA_PRS_bmi_UNadj) <- c("SUBJECT_ID", "OSA_PRS_bmi_UNadj")

dat <- left_join(dat, OSA_PRS_bmi_UNadj, by = "SUBJECT_ID")

dat$std_OSA_PRS_bmi_UNadj <- as.vector(scale(dat$OSA_PRS_bmi_UNadj))

OSA_PRS_bmi_adj <- fread("/Data/OSA_PRS/bmi_adj/OSA_bmi_adj.txt", header = T)

OSA_PRS_bmi_adj <- OSA_PRS_bmi_adj[,-c("FID")]

colnames(OSA_PRS_bmi_adj) <- c("SUBJECT_ID", "OSA_PRS_bmi_adj")

dat <- left_join(dat, OSA_PRS_bmi_adj, by = "SUBJECT_ID")

dat$std_OSA_PRS_bmi_adj <- as.vector(scale(dat$OSA_PRS_bmi_adj))

outcomes <- c("Mod_severe_OSA","Mild_severe_OSA")

PRSs <- c("std_OSA_PRS_bmi_UNadj", "std_OSA_PRS_bmi_adj")

#re-code vals from Yes/No to 1/0
dat$Mod_severe_OSA <- ifelse(dat$Mod_severe_OSA == "Yes",1, 0)
dat$Mild_severe_OSA <- ifelse(dat$Mild_severe_OSA == "Yes", 1, 0)


# for each outcome, compute the association of OSA-PRS with OSA at baseline


res <- c()

for (i in 1:length(outcomes)){
  
  
  outcome <- outcomes[i]
  
  #keep only inds with PRSs available
  dat$keep <- ifelse( !is.na(dat$PRSsum_std),1, 0)
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    
    
    cur_prs <- PRSs[j]
    cur_formula <- paste(outcome, "~ Age + Gender + BMI + Center +",  
                         "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + ", cur_prs)
    
    
    mod <- svyglm(as.formula(cur_formula),
                  design = survey_obj, 
                  family = quasibinomial(link = "logit"))
    
    
    res <- rbind(res, data.frame(
      Outcome = outcome,
      n = nrow(model.matrix(mod)),
      mean_age = svymean(~Age, survey_obj),
      PRS = cur_prs,
      log_OR = summary(mod)$coef[cur_prs,"Estimate"],
      log_OR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
      log_OR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"]))
    
    
  }
}

write.csv(res, file = "/Results/assoc_OSA_PRS_baseline_OSA.csv")


