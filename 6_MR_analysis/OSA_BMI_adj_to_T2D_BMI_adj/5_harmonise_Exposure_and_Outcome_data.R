library(TwoSampleMR)


file_path <- "/OSA_BMI_adj_to_T2D_BMI_adj/"

p_value_dir_path <- "SNPs_filtered_by_p/10_to_negative_7/"
#5_times_10_to_negative_8 
#10_to_negative_5 
#10_to_negative_7 
file_name <- "merged_Exposure_and_Outcome_GWAS_after_clump.csv"

output_file <- paste0(file_path, p_value_dir_path, "harmonized_exposure_to_outcome.csv")

merged_Exposure_and_Outcome_GWAS <- fread(paste0(file_path, p_value_dir_path, file_name))

#Exposure GWAS
#subset to only the columns we need
Exposure_GWAS <- merged_Exposure_and_Outcome_GWAS[, c("SNP", "beta.exposure",
                                                  "se.exposure", "effect_allele.exposure",
                                                  "other_allele.exposure", "eaf.exposure")]


Outcome_GWAS <- merged_Exposure_and_Outcome_GWAS[, c("SNP", "beta.outcome",
                                                     "se.outcome", "effect_allele.outcome",
                                                     "other_allele.outcome", "eaf.outcome")]

#add the required columns
Exposure_GWAS$exposure <- "OSA_BMI_adj"
Exposure_GWAS$id.exposure <- "OSA"

Outcome_GWAS$outcome <- "T2D_BMI_adj"
Outcome_GWAS$id.outcome <- "T2D"

harmonized_Exposure_and_Outcome_GWASs <- harmonise_data(Exposure_GWAS, Outcome_GWAS, action = 2)

#save harmonized data
write.csv(harmonized_Exposure_and_Outcome_GWASs, output_file, row.names = FALSE)


