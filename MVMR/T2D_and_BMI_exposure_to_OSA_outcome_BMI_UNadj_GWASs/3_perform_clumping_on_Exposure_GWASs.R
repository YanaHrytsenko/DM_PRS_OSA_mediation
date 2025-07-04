library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ieugwasr)

#original GWAS before filtering by p-value
T2DM <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_T2D_BMI_UNadj_for_extract_exposure.csv")
BMI <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/Exposure_BMI_for_extract_exposure.csv")

#path to filtered files
file_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_and_BMI_exposure_to_OSA_outcome_BMI_UNadj_GWASs/SNPs_filtered_by_p/"

clump_threshold_path <- "5_times_10_to_negative_8/"
#clump_threshold_path <- "10_to_negative_7/"
#clump_threshold_path <- "10_to_negative_5/"


threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
#threshold <- "10_to_negative_5"

file_name_exposure_1 <- paste0("Exposure_T2D_BMI_UNadj_distinct_SNPs_filtered_by_p_",threshold,".csv")
file_name_exposure_2 <- paste0("Exposure_BMI_distinct_SNPs_filtered_by_p_",threshold,".csv")

path_to_file_w_new_colnames <- paste0(file_path,clump_threshold_path, "Exposure_T2D_BMI_distinct_SNPs_filtered_by_p_",threshold,"_cols_for_TwoSampleMR.csv")

output_path <- paste0(file_path,clump_threshold_path,"Exposure_T2D_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped.csv")

sig_T2DM <- fread(paste0(file_path,clump_threshold_path,file_name_exposure_1))
sig_BMI <- fread(paste0(file_path,clump_threshold_path,file_name_exposure_2))

#Create the union of all rsIDs
all_snps <- unique(c(sig_T2DM$rsID, sig_BMI$rsID))

#save the unique rsIDs
#write.csv(all_snps, paste0(file_path,clump_threshold_path,"Unique_rsIDs_instruments_T2DM_BMI.csv"), row.names = FALSE)

#Extract these SNPs from each GWAS
T2DM_filtered <- T2DM[rsID %in% all_snps]
BMI_filtered <- BMI[rsID %in% all_snps]

#reorder by rsID to make sure they follow the same order
T2DM_filtered_reordered <- T2DM_filtered[order(T2DM_filtered$rsID), ]
BMI_filtered_reordered <- BMI_filtered[order(BMI_filtered$rsID), ]

#make sure they follow same order
identical(T2DM_filtered_reordered$rsID, BMI_filtered_reordered$rsID) #TRUE

#combine the two datasets for clumping
combined_T2DM_BMI_filtered <- rbind(T2DM_filtered_reordered, BMI_filtered_reordered)

#because we want to do clumping based on T2Dm p-values (smaller sample size than BMI and thus we want to prioritize those SNPs)
#combined_T2DM_BMI_filtered$Pvalue_T2DM <- rep(T2DM_filtered_reordered$Pvalue, times = 2) 

#rename the dataframe columns before clumping
#colnames(combined_T2DM_BMI_filtered) <- c("Phenotype", "SNP", "chr", "position", "beta", "se" ,"eaf", "effect_allele" ,"other_allele", "Pvalue_by_exposure", "pval")

#this if we did not create a dummy variable
colnames(combined_T2DM_BMI_filtered) <- c("Phenotype", "SNP", "chr", "position", "beta", "se" ,"eaf", "effect_allele" ,"other_allele", "pval")


#save into a file because I will need to read with the function which assigns the columns
write.csv(combined_T2DM_BMI_filtered, path_to_file_w_new_colnames, row.names = FALSE)

head(combined_T2DM_BMI_filtered)

#read in file w new colnames to perform clumping
Exposure_data <- read_exposure_data(path_to_file_w_new_colnames, sep = ",")
head(Exposure_data)

#perform local clumping with PLINK
clumped_Exposure_data <- ld_clump(
  dplyr::tibble(rsid = Exposure_data$SNP, pval = Exposure_data$pval.exposure, id = Exposure_data$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/PLINK_files/1kg.v3/EUR"
)



#add columns needed before saving for both exposures

clumped_Exposure_data_combined_T2DM_BMI <- combined_T2DM_BMI_filtered[SNP %in% clumped_Exposure_data$rsid]

table(clumped_Exposure_data_combined_T2DM_BMI$Phenotype)


#save a list of SNPs for the Exposure after clumping
write.csv(clumped_Exposure_data_combined_T2DM_BMI, output_path, row.names = FALSE)
