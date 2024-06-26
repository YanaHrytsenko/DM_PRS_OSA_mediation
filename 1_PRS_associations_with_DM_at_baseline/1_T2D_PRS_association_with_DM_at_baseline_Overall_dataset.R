library(ggplot2)
library(tidyverse)
library(survey)
library(tableone)


dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")

#BMI-adjusted T2D-PRS
prs_BMI_adj <- fread("/Data/EUR_BMI_adjusted_diabetes_GWAS/HCHS_BMIadj_T2D_PRS_CS.txt", header = T)

prs_BMI_adj <- prs_BMI_adj[, c("IID", "std_score_BMIadj_T2D_PRS_CS")]

colnames(prs_BMI_adj) <- c("SUBJECT_ID", "BMIadj_T2D_PRS")

dat <- left_join(dat, prs_BMI_adj, by = "SUBJECT_ID")


#comparison T2D-PRSs
PGS002308_PRSs <- read.csv("/PRS_for_comparison/PGS002308_T2D_PRS.txt", sep = '\t')

PGS002308_PRSs <- PGS002308_PRSs[,c("IID", "std_PRS")]

colnames(PGS002308_PRSs) <- c("SUBJECT_ID", "PGS002308_std_PRSs")

dat <- merge(dat, PGS002308_PRSs, by = "SUBJECT_ID")

PGS003867_PRSs <- read.csv("/PRS_for_comparison/PGS003867_T2D_PRS.txt", sep = '\t')

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

write.csv(res, file = "/Results/assoc_T2D_PRS_wwith_DM_at_baseline_overall_dataset.csv")


