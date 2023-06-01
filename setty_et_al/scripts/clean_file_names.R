# Load the Required Libraries
suppressPackageStartupMessages(library(tidyverse))

# Load unprocessed Data
raw.meta <- read.csv("../../scRNASeqData/setty_et_al/File_Rename_Metadata.tsv",
  sep = "\t"
)

# Split Over Two columns
raw.meta <- raw.meta %>% separate_rows(submitted_ftp, fastq_ftp, sep = ";")

# Remove prefix
raw.meta$fastq_ftp <- substring(raw.meta$fastq_ftp, 52)
raw.meta$submitted_ftp <- substring(raw.meta$submitted_ftp, 46)

# New File_names
raw.meta$newFileNames <- "001.fastq.gz"

# Add R1 and R2
raw.meta[raw.meta$submitted_ftp %in% grep(raw.meta$submitted_ftp, pattern = "R1", value = T), "Strand"] <- "R1"
raw.meta[is.na(raw.meta$Strand), "Strand"] <- "R2"

# Add Individual Status
raw.meta$ind <- str_sub(substring(raw.meta$library_name, 7), end = -9)

# Add Sample Info
raw.meta$sample <- paste0("S", substring(raw.meta$library_name, 16))

# Add Lane Info
raw.meta[raw.meta$ind == "P3", "Lane"] <- str_sub(substring(c(raw.meta[raw.meta$ind == "P3", ]$submitted_ftp), 34), end = -17)

# Adding cutsom Lane
raw.meta[raw.meta$ind == "P1", "Lane"] <- c("L002", "L002", "L001", "L001")
raw.meta[raw.meta$ind == "P2", "Lane"] <- c("L001", "L001", "L002", "L002", "L001", "L001", "L002", "L002", "L001", "L001", "L002", "L002")

# Create New name
raw.meta$newFileNames <- paste(raw.meta$ind, raw.meta$sample, raw.meta$Lane, raw.meta$Strand, raw.meta$newFileNames, sep = "_")


for (i in c("P1", "P2", "P3")) {
  ind.df <- raw.meta[raw.meta$ind == i, ]

  change.vec <- paste("mv", ind.df$fastq_ftp, ind.df$newFileNames, sep = " ")

  write.table(change.vec,
    paste0("../../scRNASeqData/setty_et_al/", i, "_name_change.sh"),
    col.names = F, row.names = F, quote = F
  )
}


# Clean
raw.meta <- as.data.frame(raw.meta[, c(1:7)])

# Read Manifesst file
manifest <- read.csv("../../scRNASeqData/setty_et_al/manifest_from_HCa.tsv",
  sep = "\t"
)
# Subset
manifest.sub <- manifest[, colnames(manifest) %in% c(
  "donor_organism.sex",
  "file_name", "donor_organism.organism_age"
)]

# Subset2
manifest.sub <- manifest.sub[manifest.sub$file_name %in% raw.meta$submitted_ftp, ]

# Clean
colnames(manifest.sub) <- c("submitted_ftp", "sex", "age")
manifest.sub$age <- as.integer(str_sub(manifest.sub$age, end = -6))

# Merge
prcessed.meta <- merge(raw.meta, manifest.sub, "submitted_ftp")

# Write
write.table(prcessed.meta, "../../scRNASeqData/setty_et_al/File_Rename_Metadata_Processed.tsv",
  sep = "\t", row.names = F
)
