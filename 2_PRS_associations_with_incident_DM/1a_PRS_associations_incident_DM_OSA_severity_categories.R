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




# for each strata of interest, compute the association of DM-PRS with incident DM 
# here we would limit the analysis to individuals who participated in visit 2, and use visit 2 weights. 
dat <- dat[which(!is.na(dat$WEIGHT_NORM_OVERALL_V2)),]
nrow(dat) # 8797



#define outcome variables
dat$inc_diab <- 0
dat$inc_diab[which(dat$Visit2_diabetes == "Diabetic")] <- 1
dat$inc_diab[which(is.na(dat$Visit2_diabetes == "Diabetic"))] <- NA

dat$inc_diab_pre_diab <- 0
dat$inc_diab_pre_diab[which(dat$Visit2_diabetes == "Diabetic" | dat$Visit2_diabetes == "Hyperglycemic")] <- 1
dat$inc_diab_pre_diab[which(is.na(dat$Visit2_diabetes == "Diabetic"))] <- NA



# incident diabetes among people without diabetes at baseline
res <- c()
for (i in 1:length(strata)){
  stratum <- strata[i]
  if (stratum == "All"){
    dat$keep <- ifelse( !is.na(dat$PRS_sum) & 
                          dat$DIABETES2_INDICATOR == 0, 1, 0)
  } else{
    dat$keep <- ifelse( !is.na(dat$PRS_sum) & 
                          dat[[stratum]] == "Yes"  & 
                          dat$DIABETES2_INDICATOR == 0, 1, 0)
  }
  
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_V2, data=dat)
  survey_obj <- subset(survey_obj, keep == 1)
  
  for (j in 1:length(PRSs)){
    cur_prs <- PRSs[j]
    
    mod <- svyglm(as.formula(paste("inc_diab ~ Age + Gender + BMI + Center +", 
                                   "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + ",
                                   cur_prs)),
                  offset = log(Time_between_visits), 
                  design = survey_obj, 
                  family = quasipoisson(link = "log"))
    
    
    res <- rbind(res, data.frame(
      stratum = stratum,
      n = nrow(model.matrix(mod)),
      N_case = table(mod$y)[2],
      N_control = table(mod$y)[1],
      PRS = cur_prs,
      log_IRR = summary(mod)$coef[cur_prs,"Estimate"],
      log_IRR_SE = summary(mod)$coef[cur_prs,"Std. Error"],
      log_IRR_p = summary(mod)$coef[cur_prs,"Pr(>|t|)"])
    )
    
    
    
  }
}

write.csv(res, file = "/Results/assoc_PRS_incident_DM_by_OSA_severity_categories.csv")




