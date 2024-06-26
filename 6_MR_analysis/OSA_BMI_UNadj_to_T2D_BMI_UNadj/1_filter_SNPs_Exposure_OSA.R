library(readr)
library(dplyr)
library(data.table)

OSA_BMI_UNadj_distinct_SNPs_sub <- fread("/OSA_BMI_UNadj_to_T2D_BMI_UNadj/Exposure_OSA_BMI_UNadj_distinct_SNPs_sub.csv")


head(OSA_BMI_UNadj_distinct_SNPs_sub)

#filter variants to only include the ones that pass the p-value threshold of p < 5*10^-8
# p < 5*10^-8 
# p < 10^-7 
# p < 10^-5
OSA_BMI_UNadj_distinct_SNPs_filtered_by_p <- OSA_BMI_UNadj_distinct_SNPs_sub[which(as.numeric(OSA_BMI_UNadj_distinct_SNPs_sub$PValue) < 5*10^-8),]


#save filtered SNPs of the Exposure GWAS
write.csv(OSA_BMI_UNadj_distinct_SNPs_filtered_by_p, "/SNPs_filtered_by_p/5_times_10_to_negative_8/Exposure_OSA_BMI_UNadj_distinct_SNPs_filtered_by_p_5_time_10_to_negative_8.csv", row.names = FALSE)
