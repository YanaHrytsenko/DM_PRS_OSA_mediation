library(readr)
library(dplyr)
library(data.table)

T2DM <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_T2D_BMI_UNadj_for_extract_exposure.csv")

BMI <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_BMI_for_extract_exposure.csv")


#Match SNPs between TWO Exposures to make sure that the SNPs are present in both
matched_SNPs_T2D_BMI <- intersect(T2DM$rsID, BMI$rsID)

#subset each Exposure to the SNPs that matched
T2D_GWAS_sub <- T2DM[which(T2DM$rsID %in% matched_SNPs_T2D_BMI),]
BMI_GWAS_sub <- BMI[which(BMI$rsID %in% matched_SNPs_T2D_BMI),]

#filter variants to only include the ones that pass the p-value threshold of p < 5*10^-8
# p < 5*10^-8 DONE
# p < 10^-7 
# p < 10^-5 

#threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
threshold <- "10_to_negative_5"



#select instruments from each Exposure GWAS
sig_T2DM <- T2D_GWAS_sub[which(as.numeric(T2D_GWAS_sub$Pvalue) < 10^-5),]
sig_BMI <- BMI_GWAS_sub[which(as.numeric(BMI_GWAS_sub$Pvalue) < 10^-5),]

#save filtered SNPs of the Exposure GWASs
write.csv(sig_T2DM, paste0("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/SNPs_filtered_by_p/",threshold,"/Exposure_T2D_BMI_UNadj_distinct_SNPs_filtered_by_p_",threshold,".csv"), row.names = FALSE)

write.csv(sig_BMI, paste0("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/SNPs_filtered_by_p/",threshold,"/Exposure_BMI_distinct_SNPs_filtered_by_p_",threshold,".csv"), row.names = FALSE)




