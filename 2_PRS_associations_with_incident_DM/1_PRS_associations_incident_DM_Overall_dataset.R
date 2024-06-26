library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)

dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")

prs_BMI_adj <- fread("/Data/EUR_BMI_adjusted_diabetes_GWAS/HCHS_BMIadj_T2D_PRS_CS.txt", header = T)

prs_BMI_adj <- prs_BMI_adj[, c("IID", "std_score_BMIadj_T2D_PRS_CS")]

colnames(prs_BMI_adj) <- c("SUBJECT_ID", "BMIadj_T2D_PRS")

dat <- left_join(dat, prs_BMI_adj, by = "SUBJECT_ID")

PGS002308_PRSs <- read.csv("/PRS/PRS_for_comparison/PGS002308_T2D_PRS.txt", sep = '\t')

PGS002308_PRSs <- PGS002308_PRSs[,c("IID", "std_PRS")]

colnames(PGS002308_PRSs) <- c("SUBJECT_ID", "PGS002308_std_PRSs")

dat <- merge(dat, PGS002308_PRSs, by = "SUBJECT_ID")


PGS003867_PRSs <- read.csv("/PRS/PRS_for_comparison/PGS003867_T2D_PRS.txt", sep = '\t')

PGS003867_PRSs <- PGS003867_PRSs[,c("IID", "std_PRS")]


colnames(PGS003867_PRSs) <- c("SUBJECT_ID", "PGS003867_std_PRSs")


dat <- merge(dat, PGS003867_PRSs, by = "SUBJECT_ID")

strata <- c("All")

PRSs <- c("PRSstd_sum_std",
          "PRSstd_gap_std",
          "PRSstd_mgb_std", 
          "BMIadj_T2D_PRS",
          "PGS002308_std_PRSs", 
          "PGS003867_std_PRSs")




# for each strata of interest, compute the association of DM-PRS with incident DM 
# here we would limit the analysis to individuals who participated in visit 2, and use visit 2 weights. 
dat <- dat[which(!is.na(dat$WEIGHT_NORM_OVERALL_V2)),]


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

write.csv(res, file = "/Results/assoc_PRS_incident_DM.csv")






