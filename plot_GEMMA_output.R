# Load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(stats)
library(qqman)
# install.packages("qqman")

# Specify the directory containing GEMMA output files
gemma_output_dir <- "/group/jrigrp11/juliagh/GEMMA/output/log_response_ratio_0.1geno"

# Get a list of all *assoc.txt files in the directory
gemma_output_files <- list.files(path = gemma_output_dir, pattern = "\\.assoc\\.txt$", full.names = TRUE)

# Function to perform the plotting and FDR adjustment for a single file
plot_gemma_results <- function(gemma_output_file) {
  
  # Read the GEMMA output into a data frame
  gemma_data <- read.table(gemma_output_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the trait name from the filename for use as title
  filename <- basename(gemma_output_file)
  trait_name <- str_extract(filename, "^([^_]+_[^_]+)") # Extract everything before second underscore.  Adapt the regex if your filenames are structured differently.
  
  # Calculate FDR threshold
  gemma_data <- gemma_data %>%
    mutate(FDR = p.adjust(p_wald, method = "fdr"))
  
  fdr_threshold <- 0.05  # Set your desired FDR threshold
  threshold_line <- -log10(fdr_threshold)  # This isn't strictly correct. We are plotting all points, but can draw a line at -log10(0.05)
  
  # Create the Manhattan Plot with Q vals and an FDR threshold line
  p <- ggplot(gemma_data, aes(x = ps, y = -log10(FDR))) +
    geom_point(size = 1) +
    facet_wrap(~ chr, scales = "free_x", ncol = 1) +
    labs(title = paste(trait_name, ": GEMMA Association Plot (-log10(FDR-p)"),
         x = "Genomic Position (bp)",
         y = "-log10(FDR p-value)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = threshold_line, linetype = "dashed", color = "black") # Add horizontal significance line
  
  
  # Save the plot (optional) - save it with the trait name
  plot_filename <- file.path(gemma_output_dir, paste0(trait_name, "_gemma_FDR_manhattan.pdf"))
  ggsave(plot_filename, plot = p, width = 12, height = 48, units = "in") # Save the plot
  print(paste("Plot saved to:", plot_filename))
  
}

# Loop through each file and apply the plotting function
for (gemma_output_file in gemma_output_files) {
  print(paste("Processing file:", gemma_output_file))
  plot_gemma_results(gemma_output_file)
}

## for just doing specific plots:
#plot_gemma_results("/group/jrigrp11/juliagh/GEMMA/output/GS_OX_function_01.geno_treatment_1_gemma_output.assoc.txt")
#gemma_output_files <- c("/group/jrigrp11/juliagh/GEMMA/output/DOG_function_01.geno_treatment_1_gemma_output.assoc.txt","/group/jrigrp11/juliagh/GEMMA/output/GS_OX_function_01.geno_treatment_1_gemma_output.assoc.txt")
 
## adding plots with qqman automated package
for (gemma_output_file in gemma_output_files) {
  print(paste("Processing file (qqman):", gemma_output_file))
  
  # Read the GEMMA output into a data frame
  gemma_data <- read.table(gemma_output_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Ensure column names match what qqman expects
  if (!all(c("chr", "ps", "rs", "p_wald") %in% colnames(gemma_data))) {
    print(paste("Skipping file due to missing columns:", gemma_output_file))
    next
  }
  
  # Convert chr column to numeric if necessary
  gemma_data$chr <- as.numeric(gemma_data$chr)
  
  # Extract the trait name from the filename
  filename <- basename(gemma_output_file)
  trait_name <- str_extract(filename, "^([^_]+_[^_]+)")
  
  # Generate and save the Manhattan plot
  pdf(file = file.path(gemma_output_dir, paste0(trait_name, "_gemma_qqman_manhattan.pdf")),
      width = 12, height = 6)
  
  manhattan(gemma_data, chr = "chr", bp = "ps", snp = "rs", p = "p_wald")
  
  dev.off()
  print(paste("QQman Manhattan plot saved for:", trait_name))
}
