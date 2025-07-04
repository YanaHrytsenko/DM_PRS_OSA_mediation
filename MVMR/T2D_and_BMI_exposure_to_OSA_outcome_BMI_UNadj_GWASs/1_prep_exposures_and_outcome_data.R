library(readr)
library(dplyr)
library(data.table)

#here I only select columns I will need later, rename columns across all GWASs to be the same

#T2D GWAS (Exposure)

#read in T2D BMI UNadjusted GWAS data
T2D_BMI_UNadj <- fread("/Volumes/Sofer\ Lab/Summary_statistics/2022_DIAGRAM_T2D_GWAS/DIAMANTE-EUR.sumstat.txt")


#rename columns to stay consistent
colnames(T2D_BMI_UNadj) <- c("CHR", "POS",  "SNP", "rsID",  "EA", "NEA", "EAF", "BETA", "SE", "Pvalue")


#change the format of the column so we can match on later
T2D_BMI_UNadj$SNP <- paste0(T2D_BMI_UNadj$CHR, ":", T2D_BMI_UNadj$POS) 

T2D_BMI_UNadj$Phenotype <- "T2DM"


#retain only unique SNPs
T2D_BMI_UNadj_distinct_SNPs <- T2D_BMI_UNadj %>% distinct(SNP, .keep_all=TRUE)


T2D_BMI_UNadj_distinct_SNPs_reordered_cols <- T2D_BMI_UNadj_distinct_SNPs[,c("Phenotype","rsID", "CHR", "POS", "BETA","SE","EAF","EA","NEA","Pvalue")]


#convert to Upper case the variants so it matches other GWASs file formats
T2D_BMI_UNadj_distinct_SNPs_reordered_cols$EA <- toupper(T2D_BMI_UNadj_distinct_SNPs_reordered_cols$EA)
T2D_BMI_UNadj_distinct_SNPs_reordered_cols$NEA <- toupper(T2D_BMI_UNadj_distinct_SNPs_reordered_cols$NEA)
T2D_BMI_UNadj_distinct_SNPs_reordered_cols$Pvalue <- as.numeric(T2D_BMI_UNadj_distinct_SNPs_reordered_cols$Pvalue)



#BMI GWAS (second Exposure)

BMI_GWAS <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/BMI_GWAS/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")


#rename SNP column to rsID
colnames(BMI_GWAS) <- c("CHR", "POS", "rsID", "EA", "NEA", "EAF","BETA", "SE", "Pvalue", "N")


#change the format of the column so we can match on later
BMI_GWAS$SNP <- paste0(BMI_GWAS$CHR, ":", BMI_GWAS$POS) 

BMI_GWAS$Phenotype <- "BMI"



#retain only unique SNPs
BMI_GWAS_distinct_SNPs <- BMI_GWAS %>% distinct(SNP, .keep_all=TRUE)


#reorder columns for consistency
BMI_GWAS_distinct_SNPs_reordered_cols <- BMI_GWAS_distinct_SNPs[,c("Phenotype","rsID","CHR", "POS","BETA","SE","EAF","EA","NEA","Pvalue")]

BMI_GWAS_distinct_SNPs_reordered_cols$Pvalue <- as.numeric(BMI_GWAS_distinct_SNPs_reordered_cols$Pvalue)



#read in OSA BMI UNadjusted GWAS data
OSA_BMI_UNadj <- fread("/Volumes/Sofer\ Lab/Summary_statistics/2023_MVP_OSA_GWAS_sumary_stat/UnadjBMI/2021-12-14_GWAS_OSA_UnadjustedBMI_R4_EUR_both_sex.txt.gz")


#only select columns we will use (same as in the exposure GWAS)
my_OSA_BMI_UNadj <- OSA_BMI_UNadj[, c("CHR", "rsID", "POS", "A1", "A2", "EAF", "OR", "OR_se", "PValue", "n_samples")]

#take log(OR) in case of OSA as exposure
my_OSA_BMI_UNadj$OR <- as.numeric(my_OSA_BMI_UNadj$OR)
my_OSA_BMI_UNadj$OR <- log(my_OSA_BMI_UNadj$OR)

colnames(my_OSA_BMI_UNadj) <- c("CHR", "rsID", "POS", "EA", "NEA", "EAF","BETA", "SE", "Pvalue", "N")


#add SNP column to match on with my_T2D_BMI_adj
my_OSA_BMI_UNadj$SNP <- paste0(my_OSA_BMI_UNadj$CHR, ":", my_OSA_BMI_UNadj$POS) 

#retain only unique SNPs
OSA_BMI_UNadj_distinct_SNPs <- my_OSA_BMI_UNadj %>% distinct(SNP, .keep_all=TRUE)
OSA_BMI_UNadj_distinct_SNPs$Phenotype <- "OSA"

OSA_BMI_UNadj_distinct_SNPs_reordered_cols <- OSA_BMI_UNadj_distinct_SNPs[,c("Phenotype","rsID","CHR", "POS","BETA","SE","EAF","EA","NEA","Pvalue")]

OSA_BMI_UNadj_distinct_SNPs_reordered_cols$Pvalue <- as.numeric(OSA_BMI_UNadj_distinct_SNPs_reordered_cols$Pvalue)

#save the two files
#write.csv(T2D_BMI_UNadj_distinct_SNPs_reordered_cols, "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_T2D_BMI_UNadj_for_extract_exposure.csv", row.names = FALSE)
#write.csv(BMI_GWAS_distinct_SNPs_reordered_cols, "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_BMI_for_extract_exposure.csv", row.names = FALSE)
#write.csv(OSA_BMI_UNadj_distinct_SNPs_reordered_cols, "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Outcome_OSA_for_extract_outcome.csv", row.names = FALSE)

