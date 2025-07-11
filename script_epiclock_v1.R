#set your work directory here. Make sure you copy the bed files in it before running this R script. We use a NAS storage, so the folder is 
#a network shared folder
#if you want to run our data, place in the working directory the already generated dataframes stored as TSV files. This will save lot of processing time.
setwd("//Storage/LabStorage8T/Nanopore/epiclock/R_work/escu")

#load libraries
#to open *csv, *tsv
library(readr)
#df <- read.csv("df.csv")
#to open *bed
library(data.table)
#df <- fread("df.bed", header = FALSE)
library(dplyr)
library(tidyverse)
#for liftOver
library(rtracklayer)
library(GenomicRanges)

#as some dataframes take long time to load/generate, we saved dataframes as tsv. Before generating them, we check if file exist. 
#If you want to change the input data simply delete these tsv files from the working directory.

#pileup bed results
if (file.exists("pbr.tsv")) {
  pbr <- read.delim("pbr.tsv", sep = "\t", header = TRUE)
  df <- read.delim("df.tsv", sep = "\t", header = TRUE)
} else {

print("file pbr.tsv does not exist!")
#define variables
min_coverage <- 5 #threshold for minimum reads

#list all files matching the pattern "*x.bed" where x is one or more digits
bed_files <- list.files(pattern = ".*[0-9]+.*\\.bed$")
#initialize an empty list to hold data frames
all_data <- list()

#open each bed file with methyl data unfiltered (the output from modkit)
for (file in bed_files) {
  #read the BED file (assumes no header)
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  #add the filename
  df$source_file <- file
  
  #extract the first number from the filename
  num_part <- regmatches(file, regexpr("[0-9]+", file))
  df$age <- as.numeric(num_part)
 
  # Store this data frame in the list
  all_data[[file]] <- df
}

write.table(df, "df.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 

# Combine all data into one dataframe
pbr <- do.call(rbind, all_data)

#rename columns for better visualization
names(pbr)[names(pbr) == "V1"] <- "chrom" #name of reference sequence from BAM header
names(pbr)[names(pbr) == "V2"] <- "start" #0-based start position
names(pbr)[names(pbr) == "V3"] <- "end" #0-based exclusive end position
names(pbr)[names(pbr) == "V4"] <- "modified_base_code"
names(pbr)[names(pbr) == "V5"] <- "score" #Equal to Nvalid_cov.
names(pbr)[names(pbr) == "V6"] <- "strand" #'+' for positive strand '-' for negative strand, '.' when strands are combined
names(pbr)[names(pbr) == "V7"] <- "start_c" #included for compatibility
names(pbr)[names(pbr) == "V8"] <- "end_c" #included for compatibility
names(pbr)[names(pbr) == "V9"] <- "color" #included for compatibility, always 255,0,0
names(pbr)[names(pbr) == "V10"] <- "Nvalid_cov" 
names(pbr)[names(pbr) == "V11"] <- "fraction_modified" # or ONT beta values
names(pbr)[names(pbr) == "V12"] <- "Nmod"
names(pbr)[names(pbr) == "V13"] <- "Ncanonical"
names(pbr)[names(pbr) == "V14"] <- "Nother_mod"
names(pbr)[names(pbr) == "V15"] <- "Ndelete"
names(pbr)[names(pbr) == "V16"] <- "Nfail"
names(pbr)[names(pbr) == "V17"] <- "Ndiff"
names(pbr)[names(pbr) == "V18"] <- "Nnocall"
names(pbr)[names(pbr) == "V19"] <- "sample"
names(pbr)[names(pbr) == "V20"] <- "age"
#reorder columns for better analysis
pbr <- pbr[ ,c(1, 2, 3, 4, 6, 10, 11, 5, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20)]
#keep reads with minimum coverage /reads (mincov) >=4
pbr <- pbr %>%
  filter(Nvalid_cov >= min_coverage)
#divide fraction_modified (beta values) /100 and round at 2 decimals
pbr$fraction_modified <- round(pbr$fraction_modified / 100, 2)

# Add new column by pasting column 1 and column 2 together (with _ separator)
pbr$ont_cpg <- paste(pbr[[1]], pbr[[2]], sep = "_")

# Get sorted list of all unique ages (column 2)
sorted_ages <- sort(unique(pbr[[20]]))

# Create a named vector to rename columns after pivoting
age_colnames <- as.character(sorted_ages)

write.table(pbr, "pbr.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

# Pivot and sort age columns numerically
pivot_sorted <- function(df) {
  df %>%
    select(age = 20, ont_cpg, fraction_modified) %>%
    pivot_wider(names_from = age, values_from = fraction_modified, values_fn = sum, values_fill = NA) %>%
    select(ont_cpg, all_of(as.character(sorted_ages)))  # Ensure column order
}


if (file.exists("ont_h.tsv")) {
  ont_h <- read.delim("ont_h.tsv", sep = "\t", header = TRUE)
} else {
# Filter each group
ont_h <- pbr %>% filter(modified_base_code == "h")
write.table(ont_h, "ont_h.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_m.tsv")) {
  ont_m <- read.delim("ont_m.tsv", sep = "\t", header = TRUE)
} else {
ont_m <- pbr %>% filter(modified_base_code == "m")
write.table(ont_m, "ont_m.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_hm.tsv")) {
  ont_hm <- read.delim("ont_hm.tsv", sep = "\t", header = TRUE)
} else {
hm_keys <- pbr %>%
  group_by(ont_cpg) %>%
  summarize(types = paste0(sort(unique(modified_base_code)), collapse = ""), .groups = "drop") %>%
  filter(types == "hm") %>%
  pull(ont_cpg)

ont_hm <- pbr %>% filter(ont_cpg %in% hm_keys)
write.table(ont_hm, "ont_hm.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

# Function to compute intercept and slope per row
fit_row_lm <- function(row) {
  y <- as.numeric(row[-1])
  mask <- !is.na(y)
  if (sum(mask) >= 6) {
    x <- age_values[mask]
    y <- y[mask]
    model <- lm(y ~ x)
    
    intercept <- coef(model)[1]
    slope     <- coef(model)[2]
    
    # R-squared
    r_squared <- summary(model)$r.squared
    
    # p-value of the slope
    p_value <- summary(model)$coefficients["x", "Pr(>|t|)"]
    
    # correlation coefficient
    r <- cor(x, y)
    
    c(intercept, slope, r_squared, p_value, r)
  } else {
    c(NA, NA, NA, NA, NA)
  }
}

if (file.exists("ont_h_matrix.tsv")) {
  ont_h_matrix <- read.delim("ont_h_matrix.tsv", sep = "\t", header = TRUE)
} else {
# Get sorted list of all unique ages (column 20)
sorted_ages <- sort(unique(ont_h[[20]]))
age_colnames <- as.character(sorted_ages)
ont_h_matrix  <- pivot_sorted(ont_h)
age_values <- as.numeric(colnames(ont_h_matrix)[-1])

# Apply the function row-wise (excluding ont_cpg column)
regression_stats <- t(apply(ont_h_matrix, 1, fit_row_lm))
# Add to original dataframe
ont_h_matrix$intercept <- regression_stats[, 1]
ont_h_matrix$slope <- regression_stats[, 2]
ont_h_matrix$r_squared <- regression_stats[, 3]
ont_h_matrix$p_value <- regression_stats[, 4]
ont_h_matrix$correlation <- regression_stats[, 5]
write.table(ont_h_matrix, "ont_h_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_m_matrix.tsv")) {
  ont_m_matrix <- read.delim("ont_m_matrix.tsv", sep = "\t", header = TRUE)
} else {
sorted_ages <- sort(unique(ont_m[[20]]))
age_colnames <- as.character(sorted_ages)
ont_m_matrix  <- pivot_sorted(ont_m)
age_values <- as.numeric(colnames(ont_m_matrix)[-1])
# Apply the function row-wise (excluding ont_cpg column)
regression_stats <- t(apply(ont_m_matrix, 1, fit_row_lm))
# Add to original dataframe
ont_m_matrix$intercept <- regression_stats[, 1]
ont_m_matrix$slope <- regression_stats[, 2]
ont_m_matrix$r_squared <- regression_stats[, 3]
ont_m_matrix$p_value <- regression_stats[, 4]
ont_m_matrix$correlation <- regression_stats[, 5]
write.table(ont_m_matrix, "ont_m_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_hm_matrix.tsv")) {
  ont_hm_matrix <- read.delim("ont_hm_matrix.tsv", sep = "\t", header = TRUE)
} else {
sorted_ages <- sort(unique(ont_hm[[20]]))
age_colnames <- as.character(sorted_ages)
ont_hm_matrix <- pivot_sorted(ont_hm)
age_values <- as.numeric(colnames(ont_hm_matrix)[-1])
# Apply the function row-wise (excluding ont_cpg column)
regression_stats <- t(apply(ont_hm_matrix, 1, fit_row_lm))
# Add to original dataframe
ont_hm_matrix$intercept <- regression_stats[, 1]
ont_hm_matrix$slope <- regression_stats[, 2]
ont_hm_matrix$r_squared <- regression_stats[, 3]
ont_hm_matrix$p_value <- regression_stats[, 4]
ont_hm_matrix$correlation <- regression_stats[, 5]
write.table(ont_hm_matrix, "ont_hm_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

if (file.exists("ont_h_matrix_filtered.tsv")) {
  ont_h_matrix_filtered <- read.delim("ont_h_matrix_filtered.tsv", sep = "\t", header = TRUE)
} else {
ont_h_matrix_filtered <- ont_h_matrix %>%
  filter(p_value < 0.05, r_squared >= 0.8)
write.table(ont_h_matrix_filtered, "ont_h_matrix_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_m_matrix_filtered.tsv")) {
  ont_m_matrix_filtered <- read.delim("ont_m_matrix_filtered.tsv", sep = "\t", header = TRUE)
} else {
ont_m_matrix_filtered <- ont_m_matrix %>%
  filter(p_value < 0.05, r_squared >= 0.8)
write.table(ont_m_matrix_filtered, "ont_m_matrix_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_hm_matrix_filtered.tsv")) {
  ont_hm_matrix_filtered <- read.delim("ont_hm_matrix_filtered.tsv", sep = "\t", header = TRUE)
} else {
ont_hm_matrix_filtered <- ont_hm_matrix %>%
  filter(p_value < 0.05, r_squared >= 0.8)
#export as tsv
write.table(ont_hm_matrix_filtered, "ont_hm_matrix_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

if (file.exists("ont_h_matrix_filtered_noNA.tsv")) {
  ont_h_matrix_filtered_noNA <- read.delim("ont_h_matrix_filtered_noNA.tsv", sep = "\t", header = TRUE)
} else {
  ont_h_matrix_filtered_noNA <- ont_h_matrix_filtered[complete.cases(ont_h_matrix_filtered), ]
  write.table(ont_h_matrix_filtered_noNA, "ont_h_matrix_filtered_noNA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_m_matrix_filtered_noNA.tsv")) {
  ont_m_matrix_filtered_noNA <- read.delim("ont_m_matrix_filtered_noNA.tsv", sep = "\t", header = TRUE)
} else {
  ont_m_matrix_filtered_noNA <- ont_m_matrix_filtered[complete.cases(ont_m_matrix_filtered), ]
  write.table(ont_m_matrix_filtered_noNA, "ont_m_matrix_filtered_noNA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
if (file.exists("ont_hm_matrix_filtered_noNA.tsv")) {
  ont_hm_matrix_filtered_noNA <- read.delim("ont_hm_matrix_filtered_noNA.tsv", sep = "\t", header = TRUE)
} else {
  ont_hm_matrix_filtered_noNA <- ont_hm_matrix_filtered[complete.cases(ont_hm_matrix_filtered), ]
  write.table(ont_hm_matrix_filtered_noNA, "ont_hm_matrix_filtered_noNA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

#apply ElasticNet
library(glmnet)
#this function does elastic net for files like ont_m_matrix, ont_h_matrix, ont_hm_matrix generated above.
run_elastic_net <- function(input_file, output_file, alpha = 0.5, n_columns = 7) {
   cat("Reading file:", input_file, "\n")
   # Read input TSV and keep first `n_columns` sample columns
  input_matrix <- read.delim(input_file, row.names = 1, check.names = FALSE)
  input_matrix <- input_matrix[, 1:n_columns]
 
  # Transpose: samples in rows, CpGs in columns
  beta_matrix <- t(input_matrix)
  
  # Remove columns (CpGs) with missing values
  beta_matrix <- beta_matrix[, colSums(is.na(beta_matrix)) == 0]
  
  # Extract ages from row names
  ages <- as.numeric(rownames(beta_matrix))
  
  # Check for problems
  if (any(is.na(ages))) {
    stop("Non-numeric or NA values found in sample age row names.")
  }
  stopifnot(length(ages) == nrow(beta_matrix))
  
  # Run Elastic Net with cross-validation
  cv_fit <- cv.glmnet(beta_matrix, ages, alpha = alpha)
  
  # Plot CV curve
  plot(cv_fit)
  
  # Extract coefficients at best lambda
  selected_coefs <- coef(cv_fit, s = "lambda.min")
  selected_cpgs <- rownames(selected_coefs)[which(selected_coefs != 0)]
  selected_cpgs <- setdiff(selected_cpgs, "(Intercept)")
  
  # Save results to file
  write.table(selected_cpgs, file = output_file,
              quote = FALSE, row.names = FALSE, col.names = "ont_cpg", sep = "\t")
  
  cat("Selected", length(selected_cpgs), "CpGs saved to", output_file, "\n")
  
  # Return selected CpGs (optional)
  return(selected_cpgs)
}

if (file.exists("ont_m_selected_cpgs_elasticnet.tsv")) {
  ont_m_selected_cpgs_elasticnet <- read.delim("ont_m_selected_cpgs_elasticnet.tsv", sep = "\t", header = TRUE)
  names(ont_m_selected_cpgs_elasticnet)[names(ont_m_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
} else {
run_elastic_net("ont_m_matrix.tsv", "ont_m_selected_cpgs_elasticnet.tsv")
  names(ont_m_selected_cpgs_elasticnet)[names(ont_m_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
}
if (file.exists("ont_h_selected_cpgs_elasticnet.tsv")) {
  ont_h_selected_cpgs_elasticnet <- read.delim("ont_h_selected_cpgs_elasticnet.tsv", sep = "\t", header = TRUE)
  names(ont_h_selected_cpgs_elasticnet)[names(ont_h_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
} else {
run_elastic_net("ont_h_matrix.tsv", "ont_h_selected_cpgs_elasticnet.tsv")
  names(ont_h_selected_cpgs_elasticnet)[names(ont_h_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
}
if (file.exists("ont_hm_selected_cpgs_elasticnet.tsv")) {
  ont_hm_selected_cpgs_elasticnet <- read.delim("ont_hm_selected_cpgs_elasticnet.tsv", sep = "\t", header = TRUE)
  names(ont_hm_selected_cpgs_elasticnet)[names(ont_hm_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
} else {
run_elastic_net("ont_hm_matrix.tsv", "ont_hm_selected_cpgs_elasticnet.tsv")
  names(ont_hm_selected_cpgs_elasticnet)[names(ont_hm_selected_cpgs_elasticnet) == "CpG_ID"] <- "ont_cpg"
}

#looking for these elastic net resulted CpGs in our data:
ont_m_matrix_elasticnet<- ont_m_matrix %>%
  filter(ont_cpg %in% ont_m_selected_cpgs_elasticnet$ont_cpg)
write.table(ont_m_matrix_elasticnet, "ont_m_matrix_elasticnet.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
ont_h_matrix_elasticnet<- ont_h_matrix %>%
  filter(ont_cpg %in% ont_h_selected_cpgs_elasticnet$ont_cpg)
write.table(ont_h_matrix_elasticnet, "ont_h_matrix_elasticnet.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
ont_hm_matrix_elasticnet<- ont_hm_matrix %>%
  filter(ont_cpg %in% ont_hm_selected_cpgs_elasticnet$ont_cpg)
write.table(ont_hm_matrix_elasticnet, "ont_hm_matrix_elasticnet.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#union keeping the unique CpGs  (all from Illumina) from: Horvath (353), Hannum (71), Horvath2 (391), and Weidner(102)
#for Weidner 27k platform was used => hg18
#for Hannum 450k platform was used => hg19

#liftOver for 450k
if (file.exists("lifted_450k_to_hg38.tsv")) {
  lifted_450k_to_hg38 <- read.delim("lifted_450k_to_hg38.tsv", sep = "\t", header = TRUE)
} else {
  #open illumina manifest file for 450k bedchip v1.2 https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html
  illumina_manifest_450k <- read.csv("humanmethylation450_15017482_v1-2.csv")
  #fisier total in p..a, trebuie curat
  illumina_manifest_450k <- illumina_manifest_450k[ , c(1,12,13)]
  #remove some rows
  illumina_manifest_450k <- illumina_manifest_450k[-c(1, 2, 3, 4, 5, 6, 7, 8), ]
  #rename columns
  names(illumina_manifest_450k)[names(illumina_manifest_450k) == "Illumina"] <- "IlmnID"
  names(illumina_manifest_450k)[names(illumina_manifest_450k) == "X.9"] <- "CHR"
  names(illumina_manifest_450k)[names(illumina_manifest_450k) == "X.10"] <- "MAPINFO"
  #remove rows with empty strings
  illumina_manifest_450k <- illumina_manifest_450k[!(illumina_manifest_450k$MAPINFO == "" | illumina_manifest_450k$MAPINFO == " "), ]
  #add column to match ONT file
  #make MAPINFO column as numeric
  illumina_manifest_450k$MAPINFO <- as.numeric(illumina_manifest_450k$MAPINFO)
  #add suffix chr / for liftOver make sure chromosome names start with "chr"
  illumina_manifest_450k$CHR <- paste0("chr", as.character(illumina_manifest_450k$CHR))
  #liftOver from hg18 /hg19 to hg38
  #first, create GRanges object
  gr <- GRanges(
    seqnames = illumina_manifest_450k$CHR,
    ranges = IRanges(start = illumina_manifest_450k$MAPINFO, width = 1),
    IlmnID = illumina_manifest_450k$IlmnID
  )
  #download chain file in wrkd http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/ OR hg19
  #for 450 k we need lift from hg19 to hg38 https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001567
  #load the chain file
  chain <- import.chain("hg19ToHg38.over.chain")
  #perform liftOver
  lifted_450k <- liftOver(gr, chain)
  #in case of multiple mappings filter only successfully lifted positions
  lifted_gr_450k <- unlist(lifted_450k)
  #create a df with lifted coordinates
  lifted_450k_to_hg38 <- data.frame(
    IlmnID = lifted_gr_450k$IlmnID,
    CHR_hg38 = gsub("chr", "", as.character(seqnames(lifted_gr_450k))),
    MAPINFO_hg38 = start(lifted_gr_450k),
    stringsAsFactors = FALSE
  )
  write.table(lifted_450k_to_hg38, "lifted_450k_to_hg38.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}
#liftover for 27k
if (file.exists("lifted_27k_to_hg38.tsv")) {
  lifted_27k_to_hg38 <- read.delim("lifted_27k_to_hg38.tsv", sep = "\t", header = TRUE)
} else {
  #we need to do liftOver for 27k bedchip https://support.illumina.com/downloads/humanmethylation27_product_support_files.html
  #open illumina file for 27k bedchip 
  illumina_manifest_27k <- read.csv("illumina_humanmethylation27_content.csv")
  #filter columns
  illumina_manifest_27k <- illumina_manifest_27k[ , c(1,3,4)]
  #rename columns
  names(illumina_manifest_27k)[names(illumina_manifest_27k) == "Name"] <- "IlmnID"
  names(illumina_manifest_27k)[names(illumina_manifest_27k) == "Chr"] <- "CHR"
  names(illumina_manifest_27k)[names(illumina_manifest_27k) == "MapInfo"] <- "MAPINFO"
  #remove rows with empty strings
  illumina_manifest_27k <- illumina_manifest_27k[!(illumina_manifest_27k$MAPINFO == "" | illumina_manifest_27k$MAPINFO == " "), ]
  #make MAPINFO column as numeric
  illumina_manifest_27k$MAPINFO <- as.numeric(illumina_manifest_27k$MAPINFO)
  #add suffix chr / for liftOver make sure chromosome names start with "chr"
  illumina_manifest_27k$CHR <- paste0("chr", as.character(illumina_manifest_27k$CHR))
  #liftOver from hg18 to hg38
  #first, create GRanges object
  gr <- GRanges(
    seqnames = illumina_manifest_27k$CHR,
    ranges = IRanges(start = illumina_manifest_27k$MAPINFO, width = 1),
    IlmnID = illumina_manifest_27k$IlmnID
  )
  #download http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/ OR hg19
  chain <- import.chain("hg18ToHg38.over.chain") #only for hg18 OR BUILD36 from Illumina
  #perform liftOver
  lifted_27k <- liftOver(gr, chain)
  #in case of multiple mappings filter only successfully lifted positions
  lifted_gr_27k <- unlist(lifted_27k)
  #create a df with lifted coordinates
  lifted_27k_to_hg38 <- data.frame(
    IlmnID = lifted_gr_27k$IlmnID,
    CHR_hg38 = gsub("chr", "", as.character(seqnames(lifted_gr_27k))),
    MAPINFO_hg38 = start(lifted_gr_27k),
    stringsAsFactors = FALSE
  )
  write.table(lifted_27k_to_hg38, "lifted_27k_to_hg38.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

#merge 4 clocks (all from Illumina) compatible with hg38 and ONT data
if (file.exists("c4_CpGs_hg38.tsv")) {
  c4_CpGs_hg38 <- read.delim("c4_CpGs_hg38.tsv", sep = "\t", header = TRUE)
} else {
  #open Weidner (102) data
  Weidner_102CpGs <- read_csv("Weidner_102CpGs.csv")
  #merge Weidner_102CpGs and lifted_27k_to_hg38
  Weidner102CpGs_hg38 <- merge(Weidner_102CpGs, lifted_27k_to_hg38, by = "IlmnID")
  #add column
  Weidner102CpGs_hg38$From <- "Weidner102CpGs"
  
  #liftOver for Hannum 450k = hg19
  #open Hannum (71) data
  Hannum_71CpGs <- read_csv("Hannum_71CpG.csv")
  #merge Hannum_71CpGs and lifted_450k_to_hg38
  Hannum71CpGs_hg38 <- merge(Hannum_71CpGs, lifted_450k_to_hg38, by = "IlmnID") #lifted_450k_to_hg38 is already lifted with hg19
  #add column
  Hannum71CpGs_hg38$From <- "Hannum71CpGs"
  #remove some column
  Hannum71CpGs_hg38 <- Hannum71CpGs_hg38[ , c(1,4,5,6)]
  
  #liftover for Horvath (353)
  Horvath_353CpG <- read_csv("Horvath_353CpG.csv")
  #merge Horvath_353CpG and lifted_450k_to_hg38
  Horvath353CpGs_hg38 <- merge(Horvath_353CpG, lifted_450k_to_hg38, by = "IlmnID")
  #add column
  Horvath353CpGs_hg38$From <- "Horvath353CpGs_hg19"
  #remove some column
  Horvath353CpGs_hg38 <- Horvath353CpGs_hg38[ , c(1,4,5,6)]
  
  #merge Horvath_353CpG and lifted_27k_to_hg38
  Horvath353CpGs_hg18_hg38 <- merge(Horvath_353CpG, lifted_27k_to_hg38, by = "IlmnID")
  #add column
  Horvath353CpGs_hg18_hg38$From <- "Horvath353CpGs_hg18"
  #remove some column
  Horvath353CpGs_hg18_hg38 <- Horvath353CpGs_hg18_hg38[ , c(1,4,5,6)]
  
  #liftover for Horvath (391)
  Horvath2_391CpG <- read_csv("Horvath2_391CpG.csv")
  #merge Horvath_391CpG and lifted_450k_to_hg38
  Horvath391CpGs_hg19_to_hg38 <- merge(Horvath2_391CpG, lifted_450k_to_hg38, by = "IlmnID")
  #add column
  Horvath391CpGs_hg19_to_hg38$From <- "Horvath391CpGs_hg19"
  #remove some column
  Horvath391CpGs_hg19_to_hg38 <- Horvath391CpGs_hg19_to_hg38[ , c(1,4,5,6)]
  
  #merge Horvath_391CpG and lifted_27k_to_hg38
  Horvath391CpGs_hg18_hg38 <- merge(Horvath2_391CpG, lifted_27k_to_hg38, by = "IlmnID")
  #add column
  Horvath391CpGs_hg18_hg38$From <- "Horvath391CpGs_hg18"
  #remove some column
  Horvath391CpGs_hg18_hg38 <- Horvath391CpGs_hg18_hg38[ , c(1,4,5,6)]
  
  #for Horvath2's clock we make an exception and consider MAPINFO as native for hg38
  Horvath391CpGs_hg38 <- Horvath2_391CpG
  #rename column
  names(Horvath391CpGs_hg38)[names(Horvath391CpGs_hg38) == "MAPINFO"] <- "MAPINFO_hg38"
  names(Horvath391CpGs_hg38)[names(Horvath391CpGs_hg38) == "CHR"] <- "CHR_hg38"
  #remove suffix "chr"
  Horvath391CpGs_hg38$CHR_hg38 <- gsub("chr", "", Horvath391CpGs_hg38$CHR_hg38)
  #add column
  Horvath391CpGs_hg38$From <- "Horvath391CpGs_hg38"
  
  #join vertically, use the rbind function
  c4_CpGs_hg38 <- rbind(Weidner102CpGs_hg38, Hannum71CpGs_hg38, Horvath353CpGs_hg38, Horvath353CpGs_hg18_hg38, Horvath391CpGs_hg19_to_hg38, Horvath391CpGs_hg18_hg38, Horvath391CpGs_hg38)
  #add suffix chr 
  c4_CpGs_hg38$CHR_hg38 <- paste0("chr", as.character(c4_CpGs_hg38$CHR_hg38))
  #create a column for ONT analysis "ont_cpg"
  c4_CpGs_hg38 <- c4_CpGs_hg38 %>%
    mutate(start = MAPINFO_hg38 - 1,
           ont_cpg = paste0(CHR_hg38, "_", start))
  #reorder and filter columns
  c4_CpGs_hg38 <- c4_CpGs_hg38[ , c(1,6,4)] #1770 obs /CpGs
  write.table(c4_CpGs_hg38, "c4_CpGs_hg38.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

#set union with our ONT data (reads > 5)
#for filtered data (after regression) and c4_CpGs_hg38 by ont_cpg
common_HM_filtered5_c4 <- merge(c4_CpGs_hg38, ont_hm_matrix_filtered, by = "ont_cpg") # 0 obs
common_M_filtered5_c4 <- merge(c4_CpGs_hg38, ont_m_matrix_filtered, by = "ont_cpg") # 1 obs
common_H_filtered5_c4 <- merge(c4_CpGs_hg38, ont_h_matrix_filtered, by = "ont_cpg") # 0 obs
#for UNfiltered data (before regression)
common_HM5_c4 <- merge(c4_CpGs_hg38, ont_hm_matrix, by = "ont_cpg") # 1131 obs
common_M5_c4 <- merge(c4_CpGs_hg38, ont_m_matrix, by = "ont_cpg") # 1131 obs
common_H5_c4 <- merge(c4_CpGs_hg38, ont_h_matrix, by = "ont_cpg") # 1131 obs
#for ElasticNet data
common_eHM5_c4 <- merge(c4_CpGs_hg38, ont_hm_selected_cpgs_elasticnet, by = "ont_cpg") # 0 obs
common_eM5_c4 <- merge(c4_CpGs_hg38, ont_m_selected_cpgs_elasticnet, by = "ont_cpg") # 0 obs
common_eH5_c4 <- merge(c4_CpGs_hg38, ont_h_selected_cpgs_elasticnet, by = "ont_cpg") # 0 obs



plot_methylation_heatmap <- function(ont_m_matrix_filtered, input_name = "data") {
  # Load required packages
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    install.packages("grid")
  }
  library(pheatmap)
  library(grid)
  
  # Sort by correlation (descending order)
  #sorted_data <- ont_m_matrix_filtered[order(ont_m_matrix_filtered$slope, decreasing = TRUE), ]
  #filter for NA values and by correlation desc as well
  sorted_data <- ont_m_matrix_filtered[order(ont_m_matrix_filtered$slope, decreasing = TRUE), ]
  sorted_data <- sorted_data[complete.cases(sorted_data[, c("X5", "X13", "X24", "X46", "X69", "X70", "X91")]), ]
  # Extract methylation columns
  meth_matrix <- sorted_data[, c("X5", "X13", "X24", "X46", "X69", "X70", "X91")]
  
  # Set row names to CpG identifiers (ont_cpg column)
  rownames(meth_matrix) <- sorted_data$ont_cpg
  meth_matrix <- as.matrix(meth_matrix)
  meth_matrix <- apply(meth_matrix, 2, as.numeric)
  
  # Define output filename
  file_name <- paste0("methylation_heatmap_", input_name, "_sorted_by_correlation.png")
  
  # Save to PNG with correct rendering
  png(file_name, width = 1000, height = 1200, res = 150)
  
  # Create and print the heatmap object to the PNG file
  ph <- pheatmap(meth_matrix,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 display_numbers = FALSE,
                 na_col = "white",
                 color = colorRampPalette(c("blue", "red"))(50),
                 main = paste("Methylation Heatmap -", input_name),
                 show_rownames = TRUE,
           #      fontsize_row = 6,  # Reduce font size if there are many rows
           #      fontsize_col = 10,  # Optional: Adjust column label size
           #      cellwidth = 1,    # Optional: Adjust cell width for better fit
           #      cellheight = 1,   # Optional: Adjust cell height for better fit
                 silent = TRUE)
  
  # Print the heatmap into the PNG device
  print(ph)
  
  # Turn off the device to save the image
  dev.off()
  
  message("Heatmap saved to ", file_name)
}

plot_methylation_heatmap(ont_m_matrix_filtered, input_name = "ont_m_matrix_filtered")
plot_methylation_heatmap(ont_h_matrix_filtered, input_name = "ont_h_matrix_filtered")
plot_methylation_heatmap(ont_hm_matrix_filtered, input_name = "ont_hm_matrix_filtered")

plot_methylation_lines <- function(ont_m_matrix_filtered, input_name = "data",
                                   min_slope = -Inf, min_r_squared = -Inf,
                                   min_p_value = Inf, min_correlation = -Inf) {
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
  }
  library(ggplot2)
  library(reshape2)
  
  # Apply smart slope filter
  if (min_slope >= 0) {
    slope_filter <- ont_m_matrix_filtered$slope >= min_slope
  } else {
    slope_filter <- ont_m_matrix_filtered$slope <= -min_slope
  }
  
  # Apply all filters
  filtered_data <- subset(ont_m_matrix_filtered,
                          slope_filter &
                            r_squared >= min_r_squared &
                            p_value <= min_p_value )
  
  # Filter complete cases for methylation columns
  filtered_data <- filtered_data[complete.cases(filtered_data[, c("X5", "X13", "X24", "X46", "X69", "X70", "X91")]), ]
  
  # Extract methylation columns
  meth_matrix <- filtered_data[, c("X5", "X13", "X24", "X46", "X69", "X70", "X91")]
  
  # Set row names to CpG identifiers (ont_cpg column)
  rownames(meth_matrix) <- filtered_data$ont_cpg
  meth_matrix <- as.matrix(meth_matrix)
  meth_matrix <- apply(meth_matrix, 2, as.numeric)
  
  # Melt the matrix for ggplot
  meth_df <- as.data.frame(meth_matrix)
  meth_df$ont_cpg <- rownames(meth_df)
  meth_long <- melt(meth_df, id.vars = "ont_cpg", variable.name = "Position", value.name = "Methylation")
  
  # Define output filename
  file_name <- paste0("methylation_lines_", input_name, "_filtered.png")
  
  # Save to PNG
  png(file_name, width = 1000, height = 800, res = 150)
  
  # Plot
  p <- ggplot(meth_long, aes(x = Position, y = Methylation, group = ont_cpg, color = ont_cpg)) +
    geom_line() +
    geom_point(size = 1) +
    theme_minimal() +
    theme(legend.position = "none") +  # Hide legend if too many lines
    ggtitle(paste("Methylation Lines -", input_name)) +
    xlab("Position") +
    ylab("Methylation Level")
  
  print(p)
  
  dev.off()
  
  message("Line plot saved to ", file_name)
}

plot_methylation_regression_lines <- function(ont_m_matrix_filtered, input_name = "data",
                                               min_slope = -Inf, min_r_squared = -Inf,
                                               min_p_value = Inf) {
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  
  # Apply smart slope filter
  if (min_slope >= 0) {
    slope_filter <- ont_m_matrix_filtered$slope >= min_slope
  } else {
    slope_filter <- ont_m_matrix_filtered$slope <= -min_slope
  }
  
  # Apply all filters
  filtered_data <- subset(ont_m_matrix_filtered,
                          slope_filter &
                            r_squared >= min_r_squared &
                            p_value <= min_p_value)
  
  # Define the x positions (real age values)
  positions <- c(5, 13, 24, 46, 69, 70, 91)
  
  # Build predicted values using intercept and slope
  predicted_list <- lapply(1:nrow(filtered_data), function(i) {
    intercept <- filtered_data$intercept[i]
    slope <- filtered_data$slope[i]
    cpg <- filtered_data$ont_cpg[i]
    
    data.frame(
      Position = positions,    # numeric, not factor!
      Methylation = intercept + slope * positions,
      ont_cpg = cpg
    )
  })
  
  # Combine all CpG predictions
  predicted_df <- do.call(rbind, predicted_list)
  
  # Make sure Position is numeric
  predicted_df$Position <- as.numeric(predicted_df$Position)
  
  # Define output filename
  file_name <- paste0("methylation_regression_lines_", input_name, "_filtered.png")
  
  # Save to PNG
  png(file_name, width = 2000, height = 800, res = 150)
  
  # Plot
  p <- ggplot(predicted_df, aes(x = Position, y = Methylation, group = ont_cpg, color = ont_cpg)) +
    geom_line() +
    geom_point(size = 1) +
    theme_minimal() +
    theme(legend.position = "none") +  # Hide legend
    ggtitle(paste("Methylation Regression Lines -", input_name)) +
    xlab("Age") +
    ylab("Predicted Methylation Level") +
    scale_x_continuous(breaks = positions)  # Show ticks exactly at known age points
  
  print(p)
  
  dev.off()
  
  message("Regression line plot saved to ", file_name)
}

plot_methylation_regression_lines_with_points <- function(ont_m_matrix_filtered, input_name = "data",
                                                          min_slope = 0, min_r_squared = -Inf,
                                                          min_p_value = Inf) {
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
  library(ggplot2)
  library(reshape2)
  
  # Slope filter: ascending or descending only. We separate CpGs best direct correlated with ageing from the ones inverse correlated.
  if (min_slope >= 0) {
    slope_filter <- ont_m_matrix_filtered$slope >= min_slope
  } else {
    slope_filter <- ont_m_matrix_filtered$slope <= min_slope
  }
  
  # Apply all filters
  filtered_data <- subset(ont_m_matrix_filtered,
                          slope_filter &
                            r_squared >= min_r_squared &
                            p_value <= min_p_value)
  
  # Define age positions
  positions <- c(5, 13, 24, 46, 69, 70, 91)
  
  # Get top CpGs sorted by absolute slope
  top_cpgs_df <- head(filtered_data[order(-abs(filtered_data$slope)), c("ont_cpg", "slope")], 20)
  top_cpgs <- top_cpgs_df$ont_cpg
  top_cpgs_levels <- top_cpgs[order(-abs(top_cpgs_df$slope))]  # ensure legend is ordered
  
  # Build predicted values
  predicted_list <- lapply(1:nrow(filtered_data), function(i) {
    intercept <- filtered_data$intercept[i]
    slope <- filtered_data$slope[i]
    cpg <- filtered_data$ont_cpg[i]
    legend_label <- if (cpg %in% top_cpgs) cpg else NA
    
    data.frame(
      Position = positions,
      Methylation = intercept + slope * positions,
      ont_cpg = cpg,
      legend_label = legend_label
    )
  })
  predicted_df <- do.call(rbind, predicted_list)
  
  # Factor for legend ordering
  predicted_df$legend_label <- factor(predicted_df$legend_label, levels = top_cpgs_levels)
  
  # Get actual methylation data and melt it
  actual_meth <- filtered_data[, c("5", "13", "24", "46", "69", "70", "91")]
  actual_meth$ont_cpg <- filtered_data$ont_cpg
  actual_meth$legend_label <- ifelse(actual_meth$ont_cpg %in% top_cpgs, actual_meth$ont_cpg, NA)
  actual_meth$legend_label <- factor(actual_meth$legend_label, levels = top_cpgs_levels)
  
  actual_long <- melt(actual_meth, id.vars = c("ont_cpg", "legend_label"), variable.name = "Position", value.name = "Methylation")
  actual_long$Position <- as.numeric(gsub("X", "", actual_long$Position))
  
  # Output filename
  file_name <- paste0("methylation_regression_lines_with_points_", input_name, "_filtered.png")
  
  # Save plot
  png(file_name, width = 2000, height = 800, res = 150)
  
  p <- ggplot() +
    geom_line(data = predicted_df, aes(x = Position, y = Methylation, group = ont_cpg, color = legend_label), size = 0.3) +
    geom_point(data = actual_long, aes(x = Position, y = Methylation, group = ont_cpg, color = legend_label), size = 0.5) +
    theme_minimal() +
    ggtitle(paste("Methylation Regression Lines with Actual Points -", input_name)) +
    xlab("Age") +
    ylab("Methylation Level") +
    scale_x_continuous(breaks = positions) +
    scale_color_discrete(name = "Top CpGs (sorted by slope)") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme(legend.key.size = unit(0.6, "lines"))

  print(p)
  dev.off()

  message("Plot with regression lines and actual points saved to ", file_name)
}


plot_methylation_lines(ont_m_matrix_filtered_noNA, input_name = "ont_m_matrix_filtered_noNA",
                       min_slope = 0.005,
                       min_r_squared = 0.9,
                       min_p_value = 0.01)
plot_methylation_lines(ont_h_matrix_filtered_noNA, input_name = "ont_h_matrix_filtered_noNA",
                       min_slope = 0.001,
                       min_r_squared = 0.9,
                       min_p_value = 0.01)
plot_methylation_lines(ont_hm_matrix_filtered_noNA, input_name = "ont_hm_matrix_filtered_noNA",
                       min_slope = 0.005,
                       min_r_squared = 0.9,
                       min_p_value = 0.01)


plot_methylation_regression_lines(ont_m_matrix_filtered_noNA, input_name = "ont_m_matrix_filtered_noNA",
                       min_slope = 0.003,
                       min_r_squared = 0.9,
                       min_p_value = 0.01)

plot_methylation_regression_lines(ont_h_matrix_filtered_noNA, input_name = "ont_h_matrix_filtered_noNA",
                                  min_slope = 0.003,
                                  min_r_squared = 0.9,
                                  min_p_value = 0.01)

plot_methylation_regression_lines(ont_hm_matrix_filtered_noNA, input_name = "ont_hm_matrix_filtered_noNA",
                                  min_slope = 0.003,
                                  min_r_squared = 0.9,
                                  min_p_value = 0.01)



plot_methylation_regression_lines_with_points(ont_m_matrix_filtered_noNA, input_name = "ont_m_matrix_filtered_noNa_inv",
                                  min_slope = -0.003,
                                  min_r_squared = 0.9,
                                  min_p_value = 0.01)

plot_methylation_regression_lines_with_points(ont_h_matrix_filtered_noNA, input_name = "ont_h_matrix_filtered_noNA_inv",
                                  min_slope = -0.003,
                                  min_r_squared = 0.9,
                                  min_p_value = 0.01)

plot_methylation_regression_lines_with_points(ont_hm_matrix_filtered_noNA, input_name = "ont_hm_matrix_filtered_noNA_inv",
                                  min_slope = -0.003,
                                  min_r_squared = 0.9,
                                  min_p_value = 0.01)


plot_methylation_regression_lines_with_points(ont_m_matrix_filtered_noNA, input_name = "ont_m_matrix_filtered_noNa_dir",
                                              min_slope = 0.003,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)

plot_methylation_regression_lines_with_points(ont_h_matrix_filtered_noNA, input_name = "ont_h_matrix_filtered_noNA_dir",
                                              min_slope = 0.003,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)

plot_methylation_regression_lines_with_points(ont_hm_matrix_filtered_noNA, input_name = "ont_hm_matrix_filtered_noNA_dir",
                                              min_slope = 0.003,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)


plot_methylation_regression_lines_with_points(ont_m_matrix_elasticnet, input_name = "ont_m_matrix_elasticnet",
                                              min_slope = 0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)
plot_methylation_regression_lines_with_points(ont_m_matrix_elasticnet, input_name = "ont_m_matrix_elasticnet_inv",
                                              min_slope = -0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)
plot_methylation_regression_lines_with_points(ont_h_matrix_elasticnet, input_name = "ont_h_matrix_elasticnet",
                                              min_slope = 0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)
plot_methylation_regression_lines_with_points(ont_h_matrix_elasticnet, input_name = "ont_h_matrix_elasticnet_inv",
                                              min_slope = -0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)
plot_methylation_regression_lines_with_points(ont_hm_matrix_elasticnet, input_name = "ont_hm_matrix_elasticnet",
                                              min_slope = 0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)
plot_methylation_regression_lines_with_points(ont_hm_matrix_elasticnet, input_name = "ont_hm_matrix_elasticnet_inv",
                                              min_slope = -0.001,
                                              min_r_squared = 0.9,
                                              min_p_value = 0.01)

