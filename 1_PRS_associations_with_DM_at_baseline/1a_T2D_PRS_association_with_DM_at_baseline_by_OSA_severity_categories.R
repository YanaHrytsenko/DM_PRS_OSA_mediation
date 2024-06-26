library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)


dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")


#create no OSA category
dat <- dat %>% mutate(No_OSA = factor(ifelse(Mild_OSA == "No" & Mod_severe_OSA == "No" , "Yes", "No"))) 


strata <- c("All", "No_OSA", "Mild_OSA", "Mod_severe_OSA")

PRSs <- c("PRSstd_sum_std",
          "PRSstd_gap_std",
          "PRSstd_mgb_std")




# for each strata of interest, compute the association of T2D-PRSs with DM at baseline

res <- c()
for (i in 1:length(strata)){
  stratum <- strata[i]
  if (stratum == "All"){
    dat$keep <- ifelse( !is.na(dat$PRSsum_std)  ,1, 0)
  } else{
    dat$keep <- ifelse( !is.na(dat$PRSsum_std) & 
                          dat[[stratum]] == "Yes" ,1, 0)
  }
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    cur_prs <- PRSs[j]
    cur_formula <- paste("DIABETES2_INDICATOR ~ Age + Gender + BMI + Center +",  
                         "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + ", cur_prs)
    
    
    mod <- svyglm(as.formula(cur_formula),
                  design = survey_obj, 
                  family = quasibinomial(link = "logit"))
    
    res <- rbind(res, data.frame(
      stratum = stratum,
      n = nrow(model.matrix(mod)),
      mean_age = ifelse( stratum == "All", mean(na.omit(dat$AGE)), mean(na.omit(dat$AGE[which(dat[[stratum]] == "Yes")]))),
      N_case = table(mod$y)[2],
      N_control = table(mod$y)[1],
      PRS = cur_prs,
      log_OR = summary(mod)$coef[cur_prs,"Estimate"],
      log_OR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
      log_OR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"]))
    
    
  }
}

write.csv(res, file = "/Results/assoc_T2D_PRS_wwith_DM_at_baseline_by_OSA_severity_categories.csv")


