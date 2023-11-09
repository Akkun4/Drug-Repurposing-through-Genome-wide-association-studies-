library(dplyr)
library(stringr)
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(TwoSampleMR)
library(data.table)
#opening the rds file 
#MMP_Artery_aorta_l.rds
#MMP_Artery_Appendage_l.rds
#MMP_Artery_coronary_l.rds
#MMP_Artery_tibial_l.rds
#MMP_Left_ventricle_l.rds
Open_filtered_df <-readRDS("MMP_Left_ventricle_l.rds")
# getting the unique genes
unique_genes <- unique(Open_filtered_df$hgnc_symbol) ##hgnc_symbol GENES - MMP23A, MMP7, MMP3, MMP1, MMP17 & MMP16 ****

# Counting the number of occurrences of each gene and convert to a data frame
gene_counts <- as.data.frame(table(Open_filtered_df$hgnc_symbol))

# Renaming the columns
colnames(gene_counts) <- c("Gene", "Count")

gene_counts

# Opening the GWAS files 
Gwas_formatted_files <- fread("HF_HRC_GWAS_UKBB_EUR_refGenome38.txt")


# looping through the genes to perform MR analysis
while (length(unique_genes) > 0) {
  
  # getting first gene in vector
  gene <- head(unique_genes, 1)
  
  # removing gene from vector
  unique_genes <- unique_genes[-1]
  
  # subseting data for current gene
  eQTL_MR_data <- subset(Open_filtered_df, hgnc_symbol == gene, select = c(snp_col, beta_col, se_col, eaf_col, pval_col, effect_allele_col, other_allele_col,ma_samples))
  
  # formatting data for the current gene
  MR_eqtl <- format_data(eQTL_MR_data, type = "exposure",
                         snp_col = "snp_col",
                         beta_col = "beta_col",
                         se_col = "se_col",
                         effect_allele_col = "effect_allele_col",
                         other_allele_col = "other_allele_col",
                         eaf_col='eaf_col',
                         gene_col = 'hgnc_symbol',
                         samplesize_col = "ma_samples",
                         pval_col = "pval_col") %>%
    mutate(exposure=paste(gene,"MMP_Left_ventricle"))
  
  # reading in outcome data
  outcome_dat <- read_outcome_data(
    snps = MR_eqtl$snp_col,
    filename = "HF_HRC_GWAS_UKBB_EUR_refGenome38.txt",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FRQ",
    pval_col = "P",
    units_col = "Unit",
    samplesize_col = "N"
  )
  
  # harmonizing exposure and outcome data
  dat <- harmonise_data(exposure_dat = MR_eqtl, outcome_dat = outcome_dat)
  
  # saving data to an RDS file with the gene's name and tissue 
  filename <- paste("output_files", paste(gene, "MMP_Left_ventricle", ".rds", sep = "_"), sep = "/")
  saveRDS(dat, filename)
}
# Opening harmonize RDS files. 
dat_MMP16_Artery_tibial <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP16_Artery_tibial_metallopeptidase_.rds")
datMMP16_Artery_aorta <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP16_MMP_Artery_aorta_l.rds_.rds")
datMMP28_Artery_aorta <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP28_MMP_Artery_aorta_l.rds_.rds")
datMMP11_Artery_appendages <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP11_MMP_Artery_Appendage_.rds")
datMMP16_Artery_appendages <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP16_MMP_Artery_Appendage_.rds")
datMMP16_Artery_coronarys <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP16_MMP_Artery_coronary_.rds")
datMMP11_Left_ventricle <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP11_MMP_Left_ventricle_.rds")
datMMP16_Left_ventricle <- readRDS("/Users/anrapati/Downloads/bedtoolsnew/gtex data signif_variant_gene_pairs/MMP Project/output_files/MMP16_MMP_Left_ventricle_.rds")
View(dat_MMP16_Artery_tibial)

# Heterogeneity test
mr_heterogeneity(datMMP16_Artery_aorta)
mr_heterogeneity(datMMP28_Artery_aorta)
mr_heterogeneity(dat_MMP16_Artery_tibial)
mr_heterogeneity(datMMP11_Artery_appendages)
mr_heterogeneity(datMMP16_Artery_appendages)
mr_heterogeneity(datMMP16_Artery_coronarys)
mr_heterogeneity(datMMP11_Left_ventricle)
mr_heterogeneity(datMMP16_Left_ventricle)


# Applying mr_heterogeneity() function to each dataset
result1 <- mr_heterogeneity(datMMP16_Artery_aorta)
result2 <- mr_heterogeneity(datMMP28_Artery_aorta)
result3 <- mr_heterogeneity(dat_MMP16_Artery_tibial)
result4 <- mr_heterogeneity(datMMP11_Artery_appendages)
result5 <- mr_heterogeneity(datMMP16_Artery_appendages)
result6 <- mr_heterogeneity(datMMP16_Artery_coronarys)
result7 <- mr_heterogeneity(datMMP11_Left_ventricle)
result8 <- mr_heterogeneity(datMMP16_Left_ventricle)

# Combining all the results into one table
combined_results <- rbind(result1, result2, result3, result4, result5, result6, result7, result8)
View(combined_results)
# Writing the combined results to a file
write.csv(combined_results, file = "combined_results.csv", row.names = FALSE)


#Leave one out analysis

res_loo_MMP16_Artery_aorta <- mr_leaveoneout(datMMP16_Artery_aorta)
View(res_loo_MMP16_Artery_aorta)
res_loo_MMP28_Artery_aorta <- mr_leaveoneout(datMMP28_Artery_aorta)
View(res_loo_MMP16_Artery_aorta)
res_loo_MMP16_Artery_tibial <- mr_leaveoneout(dat_MMP16_Artery_tibial)
View(res_loo_MMP16_Artery_tibial)
res_loo_MMP11_Artery_appendages <- mr_leaveoneout(datMMP11_Artery_appendages)
View(res_loo_MMP11_Artery_appendages)
res_loo_MMP16_Artery_appendages <- mr_leaveoneout(datMMP16_Artery_appendages)
View(res_loo_MMP16_Artery_appendages)
res_loo_MMP16_Artery_coronarys <- mr_leaveoneout(datMMP16_Artery_coronarys)
View(res_loo_MMP16_Artery_coronarys)
res_loo_MMP11_Left_ventricle <- mr_leaveoneout(datMMP11_Left_ventricle)
View(res_loo_MMP11_Left_ventricle)
p3 <- mr_leaveoneout_plot(res_loo_MMP11_Left_ventricle)
p3
res_loo_MMP16_Left_ventricle <- mr_leaveoneout(datMMP16_Left_ventricle)
View(res_loo_MMP16_Left_ventricle)
# Combining all the results into one table
combined_leave_out <- rbind(res_loo_MMP16_Artery_aorta, res_loo_MMP28_Artery_aorta, res_loo_MMP16_Artery_tibial,
                            res_loo_MMP11_Artery_appendages, res_loo_MMP16_Artery_appendages, res_loo_MMP16_Artery_coronarys,
                            res_loo_MMP11_Left_ventricle, res_loo_MMP16_Left_ventricle)

# Removing rows with NA values
combined_leave_out <- na.omit(combined_leave_out)
View(combined_leave_out)

# Writting the combined results to a file
write.csv(combined_results, file = "combined_leave_out.csv", row.names = FALSE)

# 2 sample MR
MR_datMMP28_Artery_aorta <- mr(datMMP28_Artery_aorta)
MR_dat_MMP16_Artery_tibial <- mr(dat_MMP16_Artery_tibial)
MR_datMMP11_Artery_appendages <- mr(datMMP11_Artery_appendages)
MR_datMMP16_Artery_appendages <- mr(datMMP16_Artery_appendages)
MR_datMMP16_Artery_coronarys <- mr(datMMP16_Artery_coronarys)
MR_datMMP11_Left_ventricle <- mr(datMMP11_Left_ventricle)
MR_datMMP16_Left_ventricle <- mr(datMMP16_Left_ventricle)
View(MR_datMMP11_Left_ventricle )

#Plots
p1  <- mr_scatter_plot(MR_datMMP11_Left_ventricle , datMMP11_Left_ventricle )
p1