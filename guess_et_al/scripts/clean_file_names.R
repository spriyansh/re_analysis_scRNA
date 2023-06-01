# Load the Required Libraries
suppressPackageStartupMessages(library(tidyverse))

# Load unprocessed Data
raw.meta <- read.csv("../../scRNASeqData/guess_et_al/ENA_Metadata.tsv",
  sep = "\t"
)
dim(raw.meta)

# Split Over Two columns
raw.meta <- raw.meta %>% separate_rows(fastq_ftp, sep = ";")

# Remove prefix
raw.meta$fastq_ftp <- substring(raw.meta$fastq_ftp, 53)

# New File_names
raw.meta$newFileNames <- "001.fastq.gz"

# Adding Patient Info
raw.meta$patient <- "P3"
raw.meta[raw.meta$sample_title %in% grep(raw.meta$sample_title, pattern = "17", value = T), "patient"] <- "P17"

# Add Disease Info
raw.meta[raw.meta$patient == "P3", "status"] <- "sAML_MDS"
raw.meta[is.na(raw.meta$status), "status"] <- c("sAML", "sAML", "MDS", "MDS")

# Demultiplexing
raw.meta[raw.meta$patient == "P3", "demultiplex"] <- "required"
raw.meta[is.na(raw.meta$demultiplex), "demultiplex"] <- "no"

# Add R1 and R2
raw.meta[raw.meta$fastq_ftp %in% grep(raw.meta$fastq_ftp, pattern = "_1.", value = T), "Strand"] <- "R1"
raw.meta[is.na(raw.meta$Strand), "Strand"] <- "R2"

# Add Sample Info
raw.meta$sample <- ifelse(raw.meta$patient == "P3", "S3", "S17")

# Add Lane Info
raw.meta$Lane <- ifelse(raw.meta$patient == "P17", "L001", "TBD")

# Create New name
raw.meta$newFileNames <- paste(paste0(raw.meta$patient, raw.meta$status), raw.meta$sample, raw.meta$Lane, raw.meta$Strand, raw.meta$newFileNames, sep = "_")

# Write file for Patient 17 only
ind.df <- raw.meta[raw.meta$patient == "P17", ]

change.vec <- paste("mv", ind.df$fastq_ftp, ind.df$newFileNames, sep = " ")
write.table(change.vec,
  paste0("../../scRNASeqData/guess_et_al/", "P17", "_name_change.sh"),
  col.names = F, row.names = F, quote = F
)


# Extract names of p13
p.13 <- raw.meta[raw.meta$patient == "P3", ]

# All Lanes
p.13$Lane <- "L001"
p.13$newFileNames <- "001.fastq.gz"

# Sample names
p.13$identify <- paste0(p.13$patient, "Comb", rep(1:8, each = 2))

# Make Filenames
p.13$newFileNames <- paste(p.13$identify, p.13$sample, p.13$Lane, p.13$Strand, p.13$newFileNames, sep = "_")


change.vec <- paste("mv", p.13$fastq_ftp, p.13$newFileNames, sep = " ")
write.table(change.vec,
  paste0("../../scRNASeqData/guess_et_al/", "P3", "_name_change.sh"),
  col.names = F, row.names = F, quote = F
)
