library(readr)
library(dplyr)
library(data.table)

T2D_BMI_adj_distinct_SNPs_sub <- fread("/T2D_BMI_adj_to_OSA_BMI_adj/Exposure_T2D_BMI_adj_distinct_SNPs_w_rsIDs.csv")


head(T2D_BMI_adj_distinct_SNPs_sub)

#filter variants to only include the ones that pass the p-value threshold of p < 5*10^-8
# p < 5*10^-8 
# p < 10^-7 
# p < 10^-5 
T2D_BMI_adj_distinct_SNPs_filtered_by_p <- T2D_BMI_adj_distinct_SNPs_sub[which(as.numeric(T2D_BMI_adj_distinct_SNPs_sub$Pvalue) < 5*10^-8),]


#save filtered SNPs of the Exposure GWAS
write.csv(T2D_BMI_adj_distinct_SNPs_filtered_by_p, "/T2D_BMI_adj_to_OSA_BMI_adj/SNPs_filtered_by_p/5_times_10_to_negative_8/Exposure_T2D_BMI_adj_distinct_SNPs_filtered_by_p_5_times_10_to_negative_8.csv", row.names = FALSE)
