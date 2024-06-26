library(TwoSampleMR)

exposure_outcome <- "OSA_BMI_adj_to_T2D_BMI_adj"

path_to_dir <- "/OSA_BMI_adj_to_T2D_BMI_adj/SNPs_filtered_by_p/"
p_value_dir_path <- "10_to_negative_5/"
#5_times_10_to_negative_8 
#10_to_negative_5 
#10_to_negative_7 #

file_name <- "harmonized_exposure_to_outcome.csv"

output_path <- "/Results/MR_analysis_results/OSA_BMI_adj_to_T2D_BMI_adj/"
output_file_1 <- paste0(output_path, p_value_dir_path, exposure_outcome, "_MR_results_TwoSampleMR.csv")

harmonized_Exposure_and_Outcome_GWASs <- fread(paste0(path_to_dir, p_value_dir_path, file_name))


# Perform MR
res <- mr(harmonized_Exposure_and_Outcome_GWASs)


write.csv(res, output_file_1, row.names = FALSE)

