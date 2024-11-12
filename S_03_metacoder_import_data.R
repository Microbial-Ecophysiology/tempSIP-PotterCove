# This script imports sequence data produced by the dada2 pipeline for use with metacoder into R

## packages ####
library(tidyverse)
library(taxa)
library(metacoder)

## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

# ASV table
asv_counts <- read.table(paste0("seq_table_for_metacoder_", prj, ".txt"))

# taxonomy table
tax_table <- read.table(paste0("tax_table_for_metacoder_", prj, ".txt"), sep = "\t")

# metadata
mdata <- read_tsv("MDATAFILE.txt")
sample.names <- mdata$Name.Mapping.File



## create taxmap object ####
obj <- parse_dada2(seq_table = asv_counts,
                   tax_table = tax_table)

# 1377 taxa (including subtaxa)
# 17,209 ASVs

## remove sample not in mapping file
obj$data$asv_table <- obj$data$asv_table[c("taxon_id", "asv_id", sample.names)]

ini_reads <- sum(obj$data$asv_table[, sample.names])
# 13,188,686 reads

## process original table ####
## remove all none 'Bacteria' sequences
obj <- filter_taxa(obj,
                   taxon_names == "Bacteria",
                   subtaxa = T)
# taxa reduced to 1320, 16,323 ASVs
sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads
# 12964394 reads -> 98.3%% of reads remain

## remove unwanted taxa Chloroplasts and Mitochondria
obj <- filter_taxa(obj, taxon_names == "Mitochondria", invert = TRUE, 
                   subtaxa = TRUE, reassign_obs = FALSE)
sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads
# 12963621 reads -> from beginning 98.3% reads remaining

obj <- filter_taxa(obj, taxon_names == "Chloroplast", invert = TRUE, 
                   subtaxa = TRUE, reassign_obs = FALSE)
sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads
# 12735859 reads -> from beginning 96.6% reads remaining

obj
# 1318  taxa (including subtaxa)
# 16,220  ASV


## remove dubletons and singletons
obj$data$asv_table <- zero_low_counts(obj, "asv_table",
                                      min_count = 3,
                                      use_total = T,
                                      other_cols = T)
# Zeroing 1426 of 16220 rows with total counts less than 3


## check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 1555 empty ASVs

# remove empty ASVs
obj <- filter_obs(obj, "asv_table",
                  ! no_reads,
                  drop_taxa = T)
obj
# ASVs reduced to 14,665, taxa 1263
sum(obj$data$asv_table[, sample.names]) 
sum(obj$data$asv_table[, sample.names])/ini_reads
# 12733021 sequences -> 96.5% reads kept from beginning

print(taxonomy_table(obj))


## calculate further tables ####
## Calculate relative abundance
# relative abundance per ASV
obj$data$rel_abd <- calc_obs_props(obj, "asv_table", other_cols = T)
# relative abundance per taxon
obj$data$rel_tax_abd <- calc_taxon_abund(obj, "rel_abd")
print(obj)

# save taxmap object for other scripts
save(obj, file = paste0("taxmap_", prj, ".RData"))


### export asv table ####
export_name <- prj

# function to export OTU table from taxmap
as.otu.table <- function(object, filename) {
  o1 <- as.data.frame(object$data$save)
  o2 <- as.data.frame(classifications(object))
  o2$taxon_id <- rownames(o2)
  o3 <- dplyr::left_join(o1, o2)
  write.table(o3, filename, sep = "\t", row.names = FALSE)
  otu_tax <<- o3
}

# name object table as object 'save' for exporting
obj$data$save <- obj$data$asv_table
as.otu.table(obj, paste0("asv_count_table_curated_", export_name, ".txt"))

obj$data$save <- obj$data$rel_abd
as.otu.table(obj, paste0("rel_abd_table_", export_name, ".txt"))

obj$data$save <- obj$data$rel_tax_abd
as.otu.table(obj, paste0("rel_tax_abd_table_", export_name, ".txt"))



# prepare table for dbRDA analysis ####
library(tidyverse)
library(taxa)
library(metacoder)

setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

# metadata
mdata <- read_tsv("MDATAFILE.txt")
sample.names <- mdata$Name.Mapping.File

load(paste0("taxmap_", prj, ".RData"))

rel.abd_forsubset <- as.data.frame(t(obj$data$rel_abd[, sample.names]))
save(rel.abd_forsubset, file = paste0(prj, "_for_dbRDA_asv-table_final.RData"))
