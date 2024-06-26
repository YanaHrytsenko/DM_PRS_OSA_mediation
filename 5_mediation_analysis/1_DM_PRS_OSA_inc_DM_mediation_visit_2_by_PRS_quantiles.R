library(mediation)
library(lme4)
library(dplyr)
library(survey)

dat <- read.csv("/Data/cleaned_DM_PRS_data_w_covars.csv")


# here we would limit the analysis to individuals who participated in visit 2, and use visit 2 weights. 
dat <- dat[which(!is.na(dat$WEIGHT_NORM_OVERALL_V2)),] 


#"PRSstd_sum_std"
#"PRSstd_gap_std" 
PRS <- "PRSstd_mgb_std"



#define outcome variables
dat$inc_diab <- 0
dat$inc_diab[which(dat$Visit2_diabetes == "Diabetic")] <- 1
dat$inc_diab <- as.numeric(dat$inc_diab)

# recode Mild_severe_OSA and Mod_severe_OSA into 0/1
dat$Mild_severe_OSA_numeric <- case_when(dat$Mild_severe_OSA == "Yes" ~ 1,
                                         dat$Mild_severe_OSA == "No" ~ 0,
                                         is.na(dat$Mild_severe_OSA) ~ NA) 

dat$Mod_severe_OSA_numeric <- case_when(dat$Mod_severe_OSA == "Yes" ~ 1,
                                        dat$Mod_severe_OSA == "No" ~ 0,
                                        is.na(dat$Mod_severe_OSA) ~ NA) 

#define log(REI) to use a continuous target variable (we classify OSA by REI index)
dat$log_REI <- log(dat$SLPA54 + 0.5)

mediators <-  c( "Mild_severe_OSA_numeric", "Mod_severe_OSA_numeric", "log_REI")


#get quantile values for PRS
PRS_q <- quantile(dat[[PRS]], p = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)

res <- c()

for (i in 1:length(mediators)){
  
  mediator <- mediators[i]
  
  cur_prs <- PRS
  
  dat$keep <- ifelse(!is.na(dat[[cur_prs]]) & 
                       dat$DIABETES2_INDICATOR == 0 & 
                       !is.na(dat[[mediator]]) & 
                       (dat$Visit2_diabetes != "Missing_visit2_diabetes"), 1, 0)
  
  
  #create a survey study object to be used in HCHS/SOL
  survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_V2, data=dat)
  
  
  # First: fit mediator-mediator model
  if (mediator != "log_REI"){
    
    model.em <- svyglm(as.formula(paste(mediator, "~ Age + Gender + ",
                                        "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +",
                                        cur_prs)),
                       design = survey_obj,
                       subset = (keep == 1),
                       family = binomial(link = "logit"))
    
  } else{ #continuous target variable case log(REI)
    
    model.em <- svyglm(as.formula(paste(mediator, "~ Age + Gender + ",
                                        "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +",
                                        cur_prs)),
                       design = survey_obj,
                       subset = (keep == 1)) 
  }
  
  
  # Second: fit mediator + mediator-outcome model 
  model.emy <- svyglm(as.formula(paste("inc_diab", "~ Age + ", mediator, " + Gender +",
                                       "PC_1 + PC_2 + PC_3 + PC_4 + PC_5 +",
                                       cur_prs)),
                      offset = Time_between_visits,
                      design = survey_obj,
                      subset = (keep == 1),
                      family = binomial(link = "logit"))
  
  #set control values (by PRS quantiles)
  for (cntrl_idx in 1:4){
    
    
    #set treat value (by PRS quantiles)
    for (trt_idx in (cntrl_idx+1):5){
      
      
      
      # Third: use the mediation R package to fit causal mediation analysis
      med.out <- mediate(model.m = model.em, model.y = model.emy, 
                         treat = cur_prs, mediator = mediator, 
                         control.value = PRS_q[cntrl_idx], treat.value = PRS_q[trt_idx],
                         robustSE = TRUE, sims = 1000)
      
      
      res <- rbind(res, data.frame(
        baseline_PRS_quantile = names(PRS_q)[cntrl_idx],
        higher_PRS_quantile = names(PRS_q)[trt_idx],
        mediator = mediator,
        n = nrow(model.matrix(model.emy)),
        N_case_mediator = table(model.em$y)[2], 
        N_control_mediator = table(model.em$y)[1],
        N_case_outcome = table(model.emy$y)[2],
        N_control_outcome = table(model.emy$y)[1],
        PRS = cur_prs,
        ACME_Estimate = med.out$d.avg,
        ACME_p_value = med.out$d.avg.p,
        ACME_lower_CI = med.out$d.avg.ci[1],
        ACME_upper_CI = med.out$d.avg.ci[2],
        ADE_Estimate = med.out$z.avg,
        ADE_p_value = med.out$z.avg.p, 
        ADE_lower_CI= med.out$z.avg.ci[1], 
        ADE_upper_CI = med.out$z.avg.ci[2],
        prop_mediated = med.out$n.avg,
        prop_mediated_p_value = med.out$n.avg.p, 
        prop_mediated_lower_CI = med.out$n.avg.ci[1], 
        prop_mediated_upper_CI = med.out$n.avg.ci[2]))
      
    }
    
  }
  
  
  
}



write.csv(res, file = paste0("/Results/DM_PRS_OSA_inc_DM_mediation_inds_visit_2_for_", PRS, "_quantiles.csv"), row.names = FALSE)
