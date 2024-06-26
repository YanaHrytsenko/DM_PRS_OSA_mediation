library(dplyr)

#define base path for easy reading
base_path <- "/OSA_BMI_UNadj_to_T2D_BMI_UNadj/"

p_value_dir <- "SNPs_filtered_by_p/10_to_negative_7/"
#5_times_10_to_negative_8/ 
#10_to_negative_5/ 
#10_to_negative_7/ 

exposure_file <- "Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_10_to_negative_7_clumped.csv"
#"Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_5_times_10_to_negative_8_clumped.csv" 
#"Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_10_to_negative_5_clumped.csv" 
#"Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_10_to_negative_7_clumped.csv" 

outcome_file <- "Outcome_T2D_BMI_UNadj_distinct_SNPs_sub.csv"

output_file <- "merged_Exposure_and_Outcome_GWAS_after_clump.csv"

output_path <- paste0(base_path, p_value_dir, output_file)


#read in list of clumped SNPs for exposure
clumped_Exposure_GWAS <- fread(paste0(base_path, p_value_dir, exposure_file))

#read in Outcome data with intersected SNPs in Exposure GWAS
Outcome_GWAS <- fread(paste0(base_path,outcome_file))


#only select columns we need
clumped_Exposure_GWAS_sub <- clumped_Exposure_GWAS[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure")]


#only select columns we need
Outcome_GWAS_sub <- Outcome_GWAS[, c("rsID", "EA", "NEA", "EAF", "Beta", "SE")]

colnames(Outcome_GWAS_sub) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome")

#merge Exposure and Outcome GWAS on SNP column
merged_Exposure_and_Outcome_GWAS <- merge(clumped_Exposure_GWAS_sub, Outcome_GWAS_sub, by = "SNP")

#save the merged file to proceed with MR 
write.csv(merged_Exposure_and_Outcome_GWAS, output_path, row.names = FALSE)

