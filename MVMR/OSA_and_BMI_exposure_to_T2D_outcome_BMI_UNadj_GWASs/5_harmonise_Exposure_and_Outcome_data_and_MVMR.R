library(TwoSampleMR)

file_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/OSA_and_BMI_exposure_to_T2D_outcome_BMI_UNadj_GWASs/"

output_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Code/Code_for_GitHub/Results/MR_analysis_results/OSA_and_BMI_exposure_to_T2D_outcome_BMI_UNadj_GWASs/"

#p_value_dir_path <- "5_times_10_to_negative_8/"
#p_value_dir_path <- "10_to_negative_7/"
p_value_dir_path <- "10_to_negative_5/"

#threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
threshold <- "10_to_negative_5"

mv_results_path <- paste0(output_path,p_value_dir_path,"MVMR_results.csv")

#read in data
Exposure <- fread(paste0(file_path,"SNPs_filtered_by_p/",p_value_dir_path,"Exposure_OSA_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped_clean_columns.csv"))
Outcome <- fread(paste0(file_path,"SNPs_filtered_by_p/",p_value_dir_path,"Outcome_GWAS_after_clump_sub_to_Exposure_SNPs.csv"))

#note: action 2 means: "Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
harmonized_Exposure_and_Outcome_GWASs <- mv_harmonise_data(Exposure, Outcome, harmonise_strictness = 2)


#5_times_10_to_negative_8:
# Harmonising OSA (OSA) and T2DM (T2DM)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10494474, rs10824218, rs10887578, rs10921744, rs10923715, rs10930656, rs10943902, rs11037497, 
# rs12584042, rs12647340, rs12987009, rs1372177, rs138289, rs1436344, rs1454687, rs1521527, rs189843, 
# rs1937433, rs2028395, rs2163188, rs2261784, rs2284746, rs2617881, rs3736328, rs6011457, rs630602, 
# rs6595205, rs6657613, rs676749, rs7129797, rs73234, rs7568228, rs7670155, rs7866532, rs8049073, rs9313052, 
# rs9929792

#10_to_negative_7
# Harmonising OSA (OSA) and T2DM (T2DM)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10824218, rs10887578, rs10921744, rs10923715, rs10930656, rs10943902, rs11037497, rs12584042, 
# rs12987009, rs1372177, rs138289, rs1436344, rs1454687, rs1521527, rs189843, rs1937433, rs2028395, 
# rs2163188, rs2284746, rs2512890, rs2617881, rs284224, rs3736328, rs5028047, rs6011457, rs630602, 
# rs6595205, rs6657613, rs676749, rs7129797, rs722290, rs73234, rs7568228, rs7670155, rs7801766, rs7866532, 
# rs8049073, rs9313052, rs9650469, rs9929792

#10_to_negative_5
# Harmonising OSA (OSA) and T2DM (T2DM)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10408139, rs10824218, rs10838439, rs10876410, rs10887578, rs11030307, rs12597712, rs12987009, 
# rs1372177, rs138289, rs1436344, rs1454687, rs1521527, rs189843, rs1937433, rs2163188, rs2267922, 
# rs2268644, rs2284746, rs2512890, rs2617881, rs284224, rs3134359, rs326202, rs379281, rs4969475, rs6011457, 
# rs6042292, rs630602, rs6595205, rs676749, rs6783231, rs7129797, rs73234, rs7568228, rs7642824, rs7866532, 
# rs9650469

mv_results <- mv_multiple(harmonized_Exposure_and_Outcome_GWASs)

write.csv(mv_results, mv_results_path, row.names = FALSE)



