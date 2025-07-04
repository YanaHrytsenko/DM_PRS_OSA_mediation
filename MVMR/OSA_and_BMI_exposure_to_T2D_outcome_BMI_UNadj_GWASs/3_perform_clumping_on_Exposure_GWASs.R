library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ieugwasr)

#original GWAS before filtering by p-value
OSA <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/T2D_BMI_UNadj_to_OSA_BMI_UNadj/Outcome_OSA_BMI_UNadj_distinct_SNPs_sub.csv")
BMI <- fread("/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/OSA_and_BMI_exposure_to_T2D_outcome_BMI_UNadj_GWASs/Exposure_BMI_for_extract_exposure.csv")


#continue
#path to filtered files
file_path <- "/Users/yana/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_diabetes_PRS_interaction_sleep/Data/post_processed_data_for_MR_analysyis/OSA_and_BMI_exposure_to_T2D_outcome_BMI_UNadj_GWASs/SNPs_filtered_by_p/"

#clump_threshold_path <- "5_times_10_to_negative_8/"
#clump_threshold_path <- "10_to_negative_7/"
clump_threshold_path <- "10_to_negative_5/"


#threshold <- "5_times_10_to_negative_8"
#threshold <- "10_to_negative_7"
threshold <- "10_to_negative_5"

file_name_exposure_1 <- paste0("Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_",threshold,".csv")
file_name_exposure_2 <- paste0("Exposure_BMI_distinct_SNPs_filtered_by_p_",threshold,".csv")

path_to_file_w_new_colnames <- paste0(file_path,clump_threshold_path, "Exposure_OSA_BMI_distinct_SNPs_filtered_by_p_",threshold,"_cols_for_TwoSampleMR.csv")

output_path <- paste0(file_path,clump_threshold_path,"Exposure_OSA_BMI_distinct_SNPs_filtered_by_p_",threshold,"_clumped.csv")

sig_OSA <- fread(paste0(file_path,clump_threshold_path,file_name_exposure_1))
sig_BMI <- fread(paste0(file_path,clump_threshold_path,file_name_exposure_2))

#Create the union of all rsIDs
all_snps <- unique(c(sig_OSA$rsID, sig_BMI$rsID))

#save the unique rsIDs
write.csv(all_snps, paste0(file_path,clump_threshold_path,"Unique_rsIDs_instruments_OSA_BMI.csv"), row.names = FALSE)

#Extract these SNPs from each GWAS
OSA_filtered <- OSA[rsID %in% all_snps]
BMI_filtered <- BMI[rsID %in% all_snps]

#reorder by rsID to make sure they follow the same order
OSA_filtered_reordered <- OSA_filtered[order(OSA_filtered$rsID), ]
BMI_filtered_reordered <- BMI_filtered[order(BMI_filtered$rsID), ]

#NOTE: the EA and NEA are flipped between BMI and OSA GWASs -- will get fixed after harmonization?
head(OSA_filtered_reordered)
head(BMI_filtered_reordered)


#make sure they follow same order
identical(OSA_filtered_reordered$rsID, BMI_filtered_reordered$rsID) #TRUE

#add phenotype column to OSA
OSA_filtered_reordered$Phenotype <- "OSA"

#take log(OR) in case of OSA as exposure
OSA_filtered_reordered$OR <- as.numeric(OSA_filtered_reordered$OR)
OSA_filtered_reordered$OR <- log(OSA_filtered_reordered$OR)

#extract only the columns we need from OSA
OSA_filtered_reordered <- OSA_filtered_reordered[,c("Phenotype","rsID", "CHR", "POS", "OR", "OR_se", "EAF", "A1", "A2", "PValue")]

#rename columns for rbind w BMI
colnames(OSA_filtered_reordered) <- c("Phenotype","rsID","CHR","POS","BETA","SE","EAF","EA","NEA","Pvalue")

#combine the two datasets for clumping
combined_OSA_BMI_filtered <- rbind(OSA_filtered_reordered, BMI_filtered_reordered)

#because we want to do clumping based on T2Dm p-values (smaller sample size than BMI and thus we want to prioritize those SNPs)
#combined_OSA_BMI_filtered$Pvalue_T2DM <- rep(OSA_filtered_reordered$PValue, times = 2) 

#rename the dataframe columns before clumping
#colnames(combined_OSA_BMI_filtered) <- c("Phenotype", "SNP", "chr", "position", "beta", "se" ,"eaf", "effect_allele" ,"other_allele", "Pvalue_by_exposure", "pval")

#this if we did not create a dummy variable
colnames(combined_OSA_BMI_filtered) <- c("Phenotype", "SNP", "chr", "position", "beta", "se" ,"eaf", "effect_allele" ,"other_allele", "pval")


#save into a file because I will need to read with the function which assigns the columns
write.csv(combined_OSA_BMI_filtered, path_to_file_w_new_colnames, row.names = FALSE)

head(combined_OSA_BMI_filtered)

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

clumped_Exposure_data_combined_OSA_BMI <- combined_OSA_BMI_filtered[SNP %in% clumped_Exposure_data$rsid]

table(clumped_Exposure_data_combined_OSA_BMI$Phenotype)


#save a list of SNPs for the Exposure after clumping
write.csv(clumped_Exposure_data_combined_OSA_BMI, output_path, row.names = FALSE)
