library(mr.raps)

exposure_outcome <- "T2D_BMI_UNadj_to_OSA_BMI_UNadj"

path_to_dir <- "/T2D_BMI_UNadj_to_OSA_BMI_UNadj/SNPs_filtered_by_p/"

p_value_dir_path <- "5_times_10_to_negative_8/"
#5_times_10_to_negative_8 
#10_to_negative_5 
#10_to_negative_7 

file_in <- "harmonized_exposure_to_outcome.csv"

harmonized_data <- fread(paste0(path_to_dir, p_value_dir_path, file_in))

harmonized_data <- as.data.frame(harmonized_data)




MR_RAPS_robust <- mr.raps.overdispersed.robust(b_exp = harmonized_data$beta.exposure,
                                               b_out = harmonized_data$beta.outcome,
                                               se_exp = harmonized_data$se.exposure,
                                               se_out = harmonized_data$se.outcome,
                                               loss.function = "tukey", 
                                               niter = 10000,
                                               suppress.warning = TRUE)


MR_RAPS_robust_res <- data.frame(Method = "MR_RAPS",
                                 Estimate = MR_RAPS_robust$beta.hat, 
                                 se = MR_RAPS_robust$beta.se,
                                 CI_L = MR_RAPS_robust$beta.hat - 1.96*MR_RAPS_robust$beta.se,
                                 CI_U = MR_RAPS_robust$beta.hat + 1.96*MR_RAPS_robust$beta.se,
                                 Pvalue = MR_RAPS_robust$beta.p.value)



output_path <- paste0("/Results/MR_analysis_results/", exposure_outcome, "/")

file_out <- paste0(output_path, p_value_dir_path, exposure_outcome, "_MR_RAPS_results.csv")

write.csv(MR_RAPS_robust_res, file = file_out, row.names = FALSE)



