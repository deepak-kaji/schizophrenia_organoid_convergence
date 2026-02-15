library(readr)
library(illuminaio)

### Read in MAnifest ###

# Path to the CSV manifest
manifest_file <- "/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/vcf/GSA-24v3-0_A2.csv"

# Read the manifest, skipping the initial header lines
manifest <- read_csv(manifest_file, skip = 7)  # skip the first 7 non-data rows

### Read in IDat Files and MAke Table with Red and Green Intensities as Columns ##$

# Path to your raw IDATs
idat_dir <- "/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/NowakowskiLabOnly_Genotyping2021QB3/raw data/"

# List green channel files
grn_files <- list.files(idat_dir, pattern="_Grn.idat$", full.names = TRUE)

# Function to read one pair of IDATs
read_pair <- function(grn_file){
  red_file <- sub("_Grn.idat$", "_Red.idat", grn_file)
  message("Reading: ", grn_file, " + ", red_file)
  
  grn <- readIDAT(grn_file)
  red <- readIDAT(red_file)
  
  # Return just the raw intensities as a list
  list(G=grn$Quants, R=red$Quants)
}

# Read all files
idat_list <- lapply(grn_files, read_pair)

### Map Intensitites to GEnomic Positions ##

# Match probe IDs between IDAT data and manifest
# We're assuming that the probe IDs in `idat_list` match those in the manifest under the column "IlmnID"

# Extract the probe IDs from the intensities (these are the rownames of the intensities)
probe_ids <- rownames(idat_list[[1]]$G)

# Now, we want to match these probe IDs with the manifest
# Ensure the probe IDs in the manifest are also rownames or in the "IlmnID" column
manifest_probes <- manifest$IlmnID

# Find the matching rows in the manifest
matching_rows <- match(probe_ids, manifest_probes)

# Filter out probes that don't have a match (for potential quality control)
valid_ids <- !is.na(matching_rows)

# Now create a new version of the manifest that matches only the valid probes
manifest_filtered <- manifest[matching_rows[valid_ids], ]

# Now combine the intensities with the manifest data (mapping probe intensities to genomic positions)
genomic_data <- data.frame(
  ProbeID = probe_ids[valid_ids],
  Chr = manifest_filtered$Chr,
  Position = manifest_filtered$MapInfo,
  AlleleA = manifest_filtered$AlleleA_ProbeSeq,
  AlleleB = manifest_filtered$AlleleB_ProbeSeq,
  Green_Intensity = idat_list[[1]]$G[valid_ids, "Mean"],   # Mean Green Intensity
  Red_Intensity = idat_list[[1]]$R[valid_ids, "Mean"],     # Mean Red Intensity
  stringsAsFactors = FALSE
)

# Check the first few rows
head(genomic_data)



