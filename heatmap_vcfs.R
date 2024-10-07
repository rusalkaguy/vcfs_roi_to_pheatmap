#!/usr/bin/env Rscript
#
# parse VCFs and create a heatmap of variants in a region
#
# Left flank example
# USAGE: Rscript heatmap_vcfs.R --list list.txt --roi chr1:161500000-161520000
#
# consider re-doing with ComplexHeatmap, which supports discrete colors directly
# https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
#
# Load necessary libraries
library(optparse)
library(vcfR)
library(dplyr)
library(tidyr)
library(tools)
library(pheatmap)

# Define options for optparse
option_list <- list(
  make_option(c("-l", "--list"), type = "character", default = "list.txt",
              help = "File containing the list of VCF files", metavar = "character"),
  make_option(c("-r", "--roi"), type = "character", 
              #default = "chr1:161500000-161520000", # left uniq flank
              #default = "chr1:161520000-161590000", # left seg
              default = "chr1:161590000-161675000", # right uniq flank
              #default = "chr1:161675000-161700000", # right uniq flank
              help = "Region of interest (format: chr:start-end)", metavar = "character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if both list file and roi were provided
if (is.null(opt$list) | is.null(opt$roi)) {
  stop("You must provide both --list and --roi arguments")
}

# Extract the genomic range (Region of Interest)
roi_info <- unlist(strsplit(opt$roi, "[:-]"))
chr_roi <- roi_info[1]
start_roi <- as.numeric(roi_info[2])
end_roi <- as.numeric(roi_info[3])

# Read list of VCF files
vcf_files <- readLines(opt$list)
vcf_files <- vcf_files[!grepl("^#", vcf_files)] # Remove comment lines

# Initialize an empty list to hold VCF data
variant_data <- list()

# Function to process each VCF file
process_vcf <- function(vcf_file) {
  vcf <- read.vcfR(vcf_file)
  
  # get file name w/o path or extensions (.g.vcf)
  # expected to be in format : SAMPLE_CONTIG_CHUNK"
  # each VCF is supposed to have ONLY 1 SAMPLE!
  vcf_sname = file_path_sans_ext(file_path_sans_ext(basename(vcf_file)))
  
  # Extract fixed fields (CHROM, POS, ID, REF, ALT, QUAL)
  fixed <- as.data.frame(getFIX(vcf))
  
  # Filter based on region of interest (ROI)
  in_roi <- fixed$CHROM == chr_roi & fixed$POS >= start_roi & fixed$POS <= end_roi
  
  # Filter based on QUAL column (removing rows with '.' in QUAL)
  qual_valid <- !is.na(fixed$QUAL)
  
  # Filter based on the presence of "<NON_REF>" in ALT
  #fixed$ALT <- gsub(",<NON_REF>", "", fixed$ALT)
  
  # Apply all filters
  filtered <- fixed[in_roi & qual_valid, ]
  
  # Extract the genotype information
  haps <- extract.haps(vcf)
  
  # Subset genotype information based on filtered variants
  filtered_haps <- haps[in_roi & qual_valid, ]
  
  # create data frame and add row(sample) info
  filtered_haps_df <- as.data.frame(t(filtered_haps))
  
  # Add sample information (row names) for easier binding
  if( ncol(haps) == 1 ) {
    # sample_contig_chunk from filename
    rownames(filtered_haps_df) <- c(vcf_sname)
  } else {
    # sample_contig_chunk from filename + rowname
    rownames(filtered_haps_df) <- paste0(vcf_sname, "_", colnames(haps))
  }
  
  return(filtered_haps_df)
}

# Process each VCF file and store results
for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  variant_data[[vcf_file]] <- process_vcf(vcf_file)
}

# Function to merge columns with ".x" and ".y" suffixes
merge_columns <- function(df) {
  # Find columns that have .x and .y suffixes
  suffix_cols <- grep("\\.x$", colnames(df), value = TRUE)
  for (col in suffix_cols) {
    # Get the base column name (without .x/.y)
    base_col <- sub("\\.x$", "", col)
    col_x <- paste0(base_col, ".x")
    col_y <- paste0(base_col, ".y")
    
    # Merge .x and .y columns using coalesce (to prefer non-NA values)
    df[[base_col]] <- coalesce(df[[col_x]], df[[col_y]])
    
    # Drop the original .x and .y columns
    df[[col_x]] <- NULL
    df[[col_y]] <- NULL
  }
  return(df)
}

# Combine genotype matrices across files by row-binding
# Align the matrices: take the union of column names and fill missing with NA
combined_matrix <- Reduce(function(x, y) {
  # Convert rownames to a column for joining
  x <- tibble::rownames_to_column(as.data.frame(x), var = "Sample")
  y <- tibble::rownames_to_column(as.data.frame(y), var = "Sample")
  
  # Perform full join on the "Sample" column
  joined <- full_join(x, y, by = "Sample")
  
  # Merge columns with .x and .y suffixes
  joined <- merge_columns(joined)
  
  # Convert the "Sample" column back to rownames
  joined <- tibble::column_to_rownames(joined, var = "Sample")
  
  return(joined)
}, variant_data)

# HACK to fix column (variant) order. Depends on all variants
# having same number of digits in CHROM and START
ordered_matrix = combined_matrix[,sort(colnames(combined_matrix))]

# Convert the matrix: set non-variant calls to NA
#combined_matrix[combined_matrix == "0/0"] <- NA

# build similarity matrix based on variant string matches
homologyDistance <- function(x) {
  distMatrix <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x)) {
      # Count number of different nucleotides between sequences i and j
      distMatrix[i, j] <- sum(x[i, ] != x[j, ], na.rm=TRUE)
    }
  }
  return(as.dist(distMatrix))  # Return as a distance matrix
}

#
# map nucleotides/variants to soemthing numeric pheatmap can handle
#
nucleotideMap <- c("REF"=0,
                   "A" = 1, "T" = 2, "G" = 3, "C" = 4, 
                   "DEL"=5, "INS"=6)

# Create a custom color palette for the nucleotides
nucleotideColors <- c("0" = "white",   # REF
                      "1" = "green",  # A
                      "2" = "red",    # T
                      "3" = "black",  # G
                      "4" = "blue",   # C
                      "5" = "gray",  # DEL
                      "6" = "orange" # INS
                      
)

numericMatrix <- apply(ordered_matrix, c(1, 2), function(x) {
  #cat(paste0("value:",x,"\n"))
  if(is.na(x)){
    #cat(paste0("\t",x,"=REF>>",nucleotideMap["REF"],"\n"))
    return(nucleotideMap["REF"])
  } else {
    if(nchar(x)>1){
      #cat(paste0("\t",x,"=INS>>",nucleotideMap["INS"],"\n"))
      return(nucleotideMap["INS"])
    } else {
      #cat(paste0("\t",x,">>",nucleotideMap[x],"\n"))
      return(nucleotideMap[x])
    }
  }
})


# now convert that ALT (character) matrix into something numeric pheatmap can handle
# Create a heatmap from the matrix
pheatmap(numericMatrix, 
         main = opt$roi,
         color=nucleotideColors,
         legend = TRUE,
         legend_breaks= nucleotideColors, 
         legend_labels = names(nucleotideMap),
         cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE
         )

