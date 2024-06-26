library(TwoSampleMR)


exposure_outcome <- "OSA_BMI_adj_to_T2D_BMI_adj"

path_to_dir <- "/OSA_BMI_adj_to_T2D_BMI_adj/SNPs_filtered_by_p/"

p_value_dir_path <- "10_to_negative_5/"
#5_times_10_to_negative_8 
#10_to_negative_5 
#10_to_negative_7

file_in <- "harmonized_exposure_to_outcome.csv"

output_path <- "/Results/MR_analysis_results/OSA_BMI_adj_to_T2D_BMI_adj/"

file_out <- paste0(output_path, p_value_dir_path, exposure_outcome, "_MR_PRESSO_results_TwoSampleMR.csv")

harmonized_data <- fread(paste0(path_to_dir, p_value_dir_path, file_in))

harmonized_data <- as.data.frame(harmonized_data)

mr_presso_results <- run_mr_presso(harmonized_data, NbDistribution = 10000, SignifThreshold = 0.05)

#save results
write.csv(mr_presso_results[[1]]$`Main MR results`, file_out, row.names = FALSE)
