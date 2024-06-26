exposure_outcome <- "T2D_BMI_UNadj_to_OSA_BMI_UNadj"

path_to_dir <- "/Results/MR_analysis_results/"

p_value_dir_path <- "5_times_10_to_negative_8/"
#5_times_10_to_negative_8 
#10_to_negative_5 
#10_to_negative_7 


file_in <- paste0(exposure_outcome,"_MR_PRESSO_results_TwoSampleMR.csv")
file_out <- paste0(exposure_outcome,"_MR_PRESSO_results_TwoSampleMR_w_CIs.csv")

results <- fread(paste0(path_to_dir, exposure_outcome, "/", p_value_dir_path, file_in))

results <- as.data.frame(results)


results$CI_low <- results$`Causal Estimate` - 1.96 * results$Sd
results$CI_high <- results$`Causal Estimate` + 1.96 * results$Sd



write.csv(results, file = paste0(path_to_dir, exposure_outcome, "/", p_value_dir_path, file_out), row.names = F)
