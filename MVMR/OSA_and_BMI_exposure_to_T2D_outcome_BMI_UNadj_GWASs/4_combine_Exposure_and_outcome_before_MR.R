library(dplyr)

#define base path for easy reading
base_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/OSA_and_BMI_exposure_to_T2D_outcome_BMI_UNadj_GWASs/"

#clump_threshold_path <- "5_times_10_to_negative_8/"
#clump_threshold_path <- "10_to_negative_7/"
clump_threshold_path <- "10_to_negative_5/"


#threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
threshold <- "10_to_negative_5"

p_value_dir <- paste0("SNPs_filtered_by_p/",clump_threshold_path)

exposure_file <- paste0("Exposure_OSA_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped.csv")
output_file_exposure <- paste0("Exposure_OSA_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped_clean_columns.csv")
output_path_for_exposure <- paste0(base_path, p_value_dir, output_file_exposure)

#this file was used in previous analyses and has rsIDs
outcome_file <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/OSA_BMI_UNadj_to_T2D_BMI_UNadj/Outcome_T2D_BMI_UNadj_distinct_SNPs_sub.csv"
output_file <- "Outcome_GWAS_after_clump_sub_to_Exposure_SNPs.csv"
output_path <- paste0(base_path, p_value_dir, output_file)

unique_SNPs_path <- paste0(base_path,"SNPs_filtered_by_p/",clump_threshold_path,"Unique_rsIDs_instruments_OSA_BMI.csv")

#read in list of clumped SNPs for exposure
clumped_Exposure_GWAS <- fread(paste0(base_path, p_value_dir, exposure_file))

clumped_Exposure_GWAS$id.exposure <- clumped_Exposure_GWAS$Phenotype

#drop erroneous p-value column and re-order

#if used dummy p-val column I named actual p-val column Pvalue_by_exposure
#clumped_Exposure_GWAS <- clumped_Exposure_GWAS[,c("SNP", "beta", "se", "effect_allele" ,"other_allele", "eaf", "Pvalue_by_exposure", "id.exposure", "Phenotype")]

clumped_Exposure_GWAS <- clumped_Exposure_GWAS[,c("SNP", "beta", "se", "effect_allele" ,"other_allele", "eaf", "pval", "id.exposure", "Phenotype")]

#rename columns for easy harmonization
colnames(clumped_Exposure_GWAS) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure" ,"other_allele.exposure", "eaf.exposure", "pval.exposure", "id.exposure", "exposure")


write.csv(clumped_Exposure_GWAS, output_path_for_exposure, row.names = FALSE)

#read in Outcome data with intersected SNPs in Exposure GWAS
Outcome_GWAS <- fread(outcome_file)


#subset to columns we need
#Outcome_GWAS_reordered <- Outcome_GWAS[,c("A1","A2","EAF", "OR_se","PValue", "rsID","beta")]
Outcome_GWAS$outcome <- "T2DM"
Outcome_GWAS$id.outcome <- "T2DM"

#reorder
Outcome_GWAS_reordered <- Outcome_GWAS[,c("rsID","Beta", "SE", "EA", "NEA", "EAF", "PValue", "id.outcome", "outcome")]

#rename columns
colnames(Outcome_GWAS_reordered) <- c("SNP", "beta.outcome","se.outcome","effect_allele.outcome", "other_allele.outcome","eaf.outcome","pval.outcome","id.outcome","outcome")

#subset outcome GWAS to the Exposure SNPs
#read in unique IDs
unique_SNPs <- fread(unique_SNPs_path)

colnames(unique_SNPs) <- "SNP"
Outcome_GWAS_sub <- Outcome_GWAS_reordered[which(Outcome_GWAS_reordered$SNP %in% unique_SNPs$SNP),]

#save subsetted outcome GWASs
write.csv(Outcome_GWAS_sub,output_path, row.names = FALSE)


