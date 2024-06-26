library(readr)
library(dplyr)
library(data.table)

OSA_BMI_adj_distinct_SNPs_sub <- fread("/OSA_BMI_adj_to_T2D_BMI_adj/Exposure_OSA_BMI_adj_distinct_SNPs_sub_w_rsIDs.csv")


head(OSA_BMI_adj_distinct_SNPs_sub)

#filter variants to only include the ones that pass the p-value threshold of p < 5*10^-8, repeat for the rest
# p < 5*10^-8 
# p < 10^-7 
# p < 10^-5
OSA_BMI_adj_distinct_SNPs_filtered_by_p <- OSA_BMI_adj_distinct_SNPs_sub[which(as.numeric(OSA_BMI_adj_distinct_SNPs_sub$PValue) < 10^-5),]


#save filtered SNPs of the Exposure GWAS
write.csv(OSA_BMI_adj_distinct_SNPs_filtered_by_p, "/OSA_BMI_adj_to_T2D_BMI_adj/SNPs_filtered_by_p/10_to_negative_5/Exposure_OSA_BMI_adj_distinct_SNPs_filtered_by_p_10_to_negative_5.csv", row.names = FALSE)
