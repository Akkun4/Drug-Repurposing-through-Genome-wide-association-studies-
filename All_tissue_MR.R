

# reading in the MMP files
MMP_files <- c(
  "MMP_Artery_aorta_l.rds",
  "MMP_Artery_Appendage_l.rds",
  "MMP_Artery_coronary_l.rds",
  "MMP_Artery_tibial_l.rds",
  "MMP_Left_ventricle_l.rds",
  "MMP_Pancreas.rds",
)


# reading in GWAS file
GWAS_file <- "/summary_data_ref38.txt"

# defining output directory
output_dir <- "output_files"

# looping through MMP files and performing MR analysis
for (MMP_file in MMP_files) {
  
  # reading in MMP file
  Open_filtered_df <- readRDS(MMP_file)
  
  # getting unique genes
  unique_genes <- unique(Open_filtered_df$gene_name)
  
  # looping through genes and performing MR analysis
  for (gene in unique_genes) {
    
    # subsetting data for current gene
    eQTL_MR_data <- subset(Open_filtered_df, gene_name == gene, select = c(snp_col, beta_col, se_col, eaf_col, pval_col, effect_allele_col, other_allele_col,ma_samples))
    
    # formatting data for current gene
    MR_eqtl <- format_data(eQTL_MR_data, type = "exposure",
                           snp_col = "snp_col",
                           beta_col = "beta_col",
                           se_col = "se_col",
                           effect_allele_col = "effect_allele_col",
                           other_allele_col = "other_allele_col",
                           eaf_col='eaf_col',
                           gene_col = 'gene_name',
                           samplesize_col = "ma_samples",
                           pval_col = "pval_col") %>%
      mutate(exposure=paste(gene, MMP_file))
    
    # reading in outcome data
    outcome_dat <- read_outcome_data(
      snps = MR_eqtl$snp_col,
      filename = GWAS_file,
      snp_col = "SNP",
      beta_col = "b",
      se_col = "se",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "freq",
      pval_col = "p",
      samplesize_col = "N"
    )
    
    # harmonizing exposure and outcome data
    dat <- harmonise_data(exposure_dat = MR_eqtl, outcome_dat = outcome_dat)
    
    # saving data to an RDS file with gene's name and tissue
    filename <- paste(output_dir, paste(gene, MMP_file, ".rds", sep = "_"), sep = "/")
    saveRDS(dat, filename)
  }
}
# Setting path to folder containing RDS files
folder_path <- "/output_files/"

# Getting a list of all RDS files in folder
rds_files <- list.files(folder_path, pattern = ".rds$")

# Looping through each RDS file
for (i in 1:length(rds_files)) {
  
  # Loading data from current RDS file
  current_file <- readRDS(paste0(folder_path, rds_files[i]))
  
  # Checking if file is empty
  if (nrow(current_file) == 0) {
    # If file is empty, print a message and move on to the next file
    cat(paste0(rds_files[i], " is empty, skipping.\n"))
    next
  }
  
  # Performing LD pruning
  LD_dat <- clump_data(current_file, clump_kb = 20, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  
  # Creating a new filename for output file
  output_filename <- paste0("LD_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Saving data as an RDS file with new filename
  saveRDS(LD_dat, file = paste0(folder_path, output_filename))
}
# Getting a list of all RDS files in the folder
rds_files <- list.files(folder_path, pattern = "^LD_.*\\.rds$")

# Looping through each RDS file
for (i in 1:length(rds_files)) {
  
  # Loading data from current RDS file
  current_file <- readRDS(paste0(folder_path, rds_files[i]))
  
  # Performing heterogeneity testing
  hetro_res <- mr_heterogeneity(current_file)
  
  # Creating a new filename for heterogeneity output file
  hetro_output_filename <- paste0("hetro_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Saving heterogeneity results as an RDS file with new filename
  saveRDS(hetro_res, file = paste0(folder_path, hetro_output_filename))
  
  # Performing MR analysis
  mr_res <- mr(current_file)
  
  # Creating a new filename for MR output file
  mr_output_filename <- paste0("mr_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Saving MR results as an RDS file with new filename
  saveRDS(mr_res, file = paste0(folder_path, mr_output_filename))
}

# Getting list of all RDS files in folder starting with "hetero_"
hetero_files <- list.files(folder_path, pattern = "^hetro.*\\.rds$")

# Initializing an empty data table to store combined hetero data
hetero_data <- data.table()

# Looping through each hetero file and combining all the data
for (file in hetero_files) {
  # Reading current file
  current_file <- readRDS(file.path(folder_path, file))
  
  # Checking if current file is empty
  if (nrow(current_file) > 0) {
    
    # Adding current file to combined hetero data
    hetero_data <- rbind(hetero_data, current_file)
  } else {
    # Printing a message indicating that current file is empty
    message(paste0("File '", file, "' is empty and will be ignored."))
  }
}
# Saving combined hetero data as an RDS file
saveRDS(hetero_data, file.path(folder_path, "combined_hetero_data.rds"))
Hetero_data <- readRDS("combined_hetero_data.rds")
View(Hetero_data)


# Getting a list of all RDS files in the folder starting with "mr_"
MR_files <- list.files(folder_path, pattern = "^mr.*\\.rds$")

# Initializing an empty data table to store combined hetero data
MR_data <- data.table()

# Looping through each hetero file and combining the data
for (file in MR_files) {
  # Reading current file
  current_file <- readRDS(file.path(folder_path, file))
  
  # Checking if current file is empty
  if (nrow(current_file) > 0) {
    
    # Adding current file to the combined hetero data
   MR_data <- rbind(MR_data, current_file)
  } else {
    # Printing a message indicating that current file is empty
    message(paste0("File '", file, "' is empty and will be ignored."))
  }
}
# Saving combined hetero data as an RDS file
saveRDS(MR_data, file.path(folder_path, "combined_MR_data.rds"))
MR_data <- readRDS("combined_MR_data.rds")
View(MR_data)

# Looping through all MR files and plotting scatter plot for corresponding LD file
for (mr_file in MR_files) {
  mmp_name <- gsub("^mr_LD_|\\.rds", "", mr_file)
  ld_file <- paste0("LD_", mmp_name, ".rds")
  
  # Checking if LD file exists and is not empty
  if (file.exists(file.path(folder_path, ld_file)) & file.info(file.path(folder_path, ld_file))$size > 0) {
    current_mr_file <- readRDS(file.path(folder_path, mr_file))
    current_ld_file <- readRDS(file.path(folder_path, ld_file))
    
    # Generating scatter plot and save as PNG
    p1 <- mr_scatter_plot(current_mr_file, current_ld_file)
    fig_name <- paste0("scatter_plot_", mmp_name, ".png")
    
    png(file.path(output_folder, fig_name))
    print(p1)
    dev.off()
  } else {
    cat("Skipping file", mr_file, "as corresponding LD file does not exist or is empty.\n")
  }
}



