library(TwoSampleMR)

#path to Exposure GWAS
file_path <- "/OSA_BMI_adj_to_T2D_BMI_adj/SNPs_filtered_by_p/"


#"10_to_negative_7/" 
#"5_times_10_to_negative_8/"
clump_threshold_path <- "10_to_negative_5/"

#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_7.csv"
#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_5_times_10_to_negative_8.csv" 
file_name <- "Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_5.csv"

#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_7_cols_for_TwoSampleMR.csv" 
#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_5_times_10_to_negative_8_cols_for_TwoSampleMR.csv" 
path_to_file_w_new_colnames <- paste0(file_path,clump_threshold_path, "Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_5_cols_for_TwoSampleMR.csv")


#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_7_clumped.csv" 
#"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_5_times_10_to_negative_8_clumped.csv" 
output_path <- paste0(file_path,clump_threshold_path,"Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_5_clumped.csv")


Exposure_GWAS <- fread(paste0(file_path, clump_threshold_path, file_name))


Exposure_GWAS <- Exposure_GWAS[,c("rsID", "CHR", "POS", "A1", "A2", "EAF", "OR", "OR_se", "PValue")]


#take log(OR) in case of OSA as exposure
Exposure_GWAS$OR <- log(Exposure_GWAS$OR)


#rename columns for further analysis
colnames(Exposure_GWAS) <- c("SNP", "chr", "position", "effect_allele", "other_allele","eaf", "beta", "se", "pval")


write.csv(Exposure_GWAS, path_to_file_w_new_colnames, row.names = FALSE)


Exposure_data <- read_exposure_data(path_to_file_w_new_colnames, sep = ",")


clumped_Exposure_data <- clump_data(Exposure_data)

#save a list of SNPs for the Exposure after clumping
write.csv(clumped_Exposure_data, output_path, row.names = FALSE)

