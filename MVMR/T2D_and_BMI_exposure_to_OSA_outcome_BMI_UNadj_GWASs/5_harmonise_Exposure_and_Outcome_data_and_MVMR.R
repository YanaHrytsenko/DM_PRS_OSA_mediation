library(TwoSampleMR)
library(data.table)
file_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/"

output_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Code/Code_for_GitHub/Results/MR_analysis_results/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/"

p_value_dir_path <- "5_times_10_to_negative_8/"
#p_value_dir_path <- "10_to_negative_7/"
#p_value_dir_path <- "10_to_negative_5/"

threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
#threshold <- "10_to_negative_5"

mv_results_path <- paste0(output_path,p_value_dir_path,"MVMR_results.csv")

#read in data
Exposure <- fread(paste0(file_path,"SNPs_filtered_by_p/",p_value_dir_path,"Exposure_T2D_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped_clean_columns.csv"))
Outcome <- fread(paste0(file_path,"SNPs_filtered_by_p/",p_value_dir_path,"Outcome_GWAS_after_clump_sub_to_Exposure_SNPs.csv"))

#note: action 2 means: "Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
harmonized_Exposure_and_Outcome_GWASs <- mv_harmonise_data(Exposure, Outcome, harmonise_strictness = 2)


#5_times_10_to_negative_8:
#w fake p-value
# Harmonising T2DM (T2DM) and OSA (OSA)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10930656, rs11753875, rs13063307, rs1796330, rs1999536, rs2028395, rs2192159, rs34341, rs4714935, rs6657613, rs703972, rs765861, rs8037894

#with "real" p-value
# Harmonising T2DM (T2DM) and OSA (OSA)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10824218, rs10887578, rs10930656, rs11753875, rs12987009, rs13063307, rs1372177, rs138289, 
# rs1436344, rs1454687, rs1521527, rs1796330, rs189843, rs1937433, rs1999536, rs2028395, rs2163188, 
# rs2192159, rs2284746, rs34341, rs4714935, rs6011457, rs6595205, rs6657613, rs676749, rs703972, rs7568228, 
# rs765861, rs8037894

#10_to_negative_7
# Harmonising T2DM (T2DM) and OSA (OSA)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs1029176, rs10824218, rs10887578, rs10930656, rs11753875, rs12987009, rs13063307, rs1363558, 
# rs1372177, rs138289, rs1436344, rs1454687, rs1521527, rs1796330, rs189843, rs1937433, rs1999536, 
# rs2163188, rs2284746, rs284224, rs2962446, rs34341, rs4714935, rs570673, rs6011457, rs6595205, rs6657613, 
# rs676749, rs703972, rs7568228, rs8037894, rs9650469

#10_to_negative_5
# Harmonising T2DM (T2DM) and OSA (OSA)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10824218, rs10842356, rs10876411, rs10887578, rs12666980, rs12987009, rs1372177, rs138289, 
# rs1436344, rs1454687, rs1521527, rs1714507, rs17537593, rs1796330, rs182705, rs189843, rs1937433, 
# rs1999536, rs2148186, rs2163188, rs2284746, rs2481750, rs284224, rs34341, rs4969475, rs566123, rs6011457, 
# rs6072085, rs6117256, rs6570526, rs6595205, rs6657613, rs676749, rs703972, rs7568228, rs7606664, rs7642824, 
# rs7726296, rs8037894, rs8046424, rs9650469

mv_results <- mv_multiple(harmonized_Exposure_and_Outcome_GWASs)

#conclusion: using "fake" p-value column (the one for the GWAS with the smaller sample size) did not have much difference on the results

#write.csv(mv_results, mv_results_path, row.names = FALSE)



