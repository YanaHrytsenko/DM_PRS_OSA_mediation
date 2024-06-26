library(readr)
library(dplyr)
library(data.table)

T2D_BMI_UNadj_distinct_SNPs_sub <- fread("/T2D_BMI_UNadj_to_OSA_BMI_UNadj/Exposure_T2D_BMI_UNadj_distinct_SNPs_sub.csv")

#filter variants to only include the ones that pass the p-value threshold of p < 5*10^-8
# p < 5*10^-8 
# p < 10^-7 
# p < 10^-5 
T2D_BMI_UNadj_distinct_SNPs_filtered_by_p <- T2D_BMI_UNadj_distinct_SNPs_sub[which(as.numeric(T2D_BMI_UNadj_distinct_SNPs_sub$PValue) < 10^-5),]


#save filtered SNPs of the Exposure GWAS
write.csv(T2D_BMI_UNadj_distinct_SNPs_filtered_by_p, "/T2D_BMI_UNadj_to_OSA_BMI_UNadj/SNPs_filtered_by_p/10_to_negative_5/Exposure_T2D_BMI_UNadj_distinct_SNPs_filtered_by_p_10_to_negative_5.csv", row.names = FALSE)
