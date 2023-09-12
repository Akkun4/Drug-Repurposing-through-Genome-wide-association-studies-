

# read in the MMP files
MMP_files <- c(
  "MMP_Artery_aorta_l.rds",
  "MMP_Artery_Appendage_l.rds",
  "MMP_Artery_coronary_l.rds",
  "MMP_Artery_tibial_l.rds",
  "MMP_Left_ventricle_l.rds",
  "MMP_Pancreas.rds",
)


# read in the GWAS file
GWAS_file <- "/summary_data_ref38.txt"

# define the output directory
output_dir <- "output_files"

# loop through the MMP files and perform MR analysis
for (MMP_file in MMP_files) {
  
  # read in the MMP file
  Open_filtered_df <- readRDS(MMP_file)
  
  # get the unique genes
  unique_genes <- unique(Open_filtered_df$gene_name)
  
  # loop through the genes and perform MR analysis
  for (gene in unique_genes) {
    
    # subset the data for the current gene
    eQTL_MR_data <- subset(Open_filtered_df, gene_name == gene, select = c(snp_col, beta_col, se_col, eaf_col, pval_col, effect_allele_col, other_allele_col,ma_samples))
    
    # format the data for the current gene
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
    
    # read in the outcome data
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
    
    # harmonize the exposure and outcome data
    dat <- harmonise_data(exposure_dat = MR_eqtl, outcome_dat = outcome_dat)
    
    # save the data to an RDS file with the gene's name and tissue
    filename <- paste(output_dir, paste(gene, MMP_file, ".rds", sep = "_"), sep = "/")
    saveRDS(dat, filename)
  }
}
# Set the path to the folder containing the RDS files
folder_path <- "/output_files/"

# Get a list of all the RDS files in the folder
rds_files <- list.files(folder_path, pattern = ".rds$")

# Loop through each RDS file
for (i in 1:length(rds_files)) {
  
  # Load the data from the current RDS file
  current_file <- readRDS(paste0(folder_path, rds_files[i]))
  
  # Check if the file is empty
  if (nrow(current_file) == 0) {
    # If the file is empty, print a message and move on to the next file
    cat(paste0(rds_files[i], " is empty, skipping.\n"))
    next
  }
  
  # Perform LD pruning
  LD_dat <- clump_data(current_file, clump_kb = 20, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  
  # Create a new filename for the output file
  output_filename <- paste0("LD_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Save the data as an RDS file with the new filename
  saveRDS(LD_dat, file = paste0(folder_path, output_filename))
}
# Get a list of all the RDS files in the folder
rds_files <- list.files(folder_path, pattern = "^LD_.*\\.rds$")

# Loop through each RDS file
for (i in 1:length(rds_files)) {
  
  # Load the data from the current RDS file
  current_file <- readRDS(paste0(folder_path, rds_files[i]))
  
  # Perform heterogeneity testing
  hetro_res <- mr_heterogeneity(current_file)
  
  # Create a new filename for the heterogeneity output file
  hetro_output_filename <- paste0("hetro_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Save the heterogeneity results as an RDS file with the new filename
  saveRDS(hetro_res, file = paste0(folder_path, hetro_output_filename))
  
  # Perform MR analysis
  mr_res <- mr(current_file)
  
  # Create a new filename for the MR output file
  mr_output_filename <- paste0("mr_", gsub(".rds", "", rds_files[i]), ".rds")
  
  # Save the MR results as an RDS file with the new filename
  saveRDS(mr_res, file = paste0(folder_path, mr_output_filename))
}

# Get a list of all the RDS files in the folder starting with "hetero_"
hetero_files <- list.files(folder_path, pattern = "^hetro.*\\.rds$")

# Initialize an empty data table to store the combined hetero data
hetero_data <- data.table()

# Loop through each hetero file and combine the data
for (file in hetero_files) {
  # Read the current file
  current_file <- readRDS(file.path(folder_path, file))
  
  # Check if the current file is empty
  if (nrow(current_file) > 0) {
    
    # Add the current file to the combined hetero data
    hetero_data <- rbind(hetero_data, current_file)
  } else {
    # Print a message indicating that the current file is empty
    message(paste0("File '", file, "' is empty and will be ignored."))
  }
}
# Save the combined hetero data as an RDS file
saveRDS(hetero_data, file.path(folder_path, "combined_hetero_data.rds"))
Hetero_data <- readRDS("combined_hetero_data.rds")
View(Hetero_data)


# Get a list of all the RDS files in the folder starting with "mr_"
MR_files <- list.files(folder_path, pattern = "^mr.*\\.rds$")

# Initialize an empty data table to store the combined hetero data
MR_data <- data.table()

# Loop through each hetero file and combine the data
for (file in MR_files) {
  # Read the current file
  current_file <- readRDS(file.path(folder_path, file))
  
  # Check if the current file is empty
  if (nrow(current_file) > 0) {
    
    # Add the current file to the combined hetero data
   MR_data <- rbind(MR_data, current_file)
  } else {
    # Print a message indicating that the current file is empty
    message(paste0("File '", file, "' is empty and will be ignored."))
  }
}
# Save the combined hetero data as an RDS file
saveRDS(MR_data, file.path(folder_path, "combined_MR_data.rds"))
MR_data <- readRDS("combined_MR_data.rds")
View(MR_data)

# Loop through all MR files and plot scatter plot for corresponding LD file
for (mr_file in MR_files) {
  mmp_name <- gsub("^mr_LD_|\\.rds", "", mr_file)
  ld_file <- paste0("LD_", mmp_name, ".rds")
  
  # Check if LD file exists and is not empty
  if (file.exists(file.path(folder_path, ld_file)) & file.info(file.path(folder_path, ld_file))$size > 0) {
    current_mr_file <- readRDS(file.path(folder_path, mr_file))
    current_ld_file <- readRDS(file.path(folder_path, ld_file))
    
    # Generate scatter plot and save as PNG
    p1 <- mr_scatter_plot(current_mr_file, current_ld_file)
    fig_name <- paste0("scatter_plot_", mmp_name, ".png")
    
    png(file.path(output_folder, fig_name))
    print(p1)
    dev.off()
  } else {
    cat("Skipping file", mr_file, "as corresponding LD file does not exist or is empty.\n")
  }
}



