# Set up libraries
if (!require("rats")) devtools::install_github("bartongroup/rats", ref="master")
library(rats)
library(tidyverse)
library(wasabi)

# Imort parameters from snakemake
cutoff <- as.numeric(snakemake@params[['cutoff']])
padj <- as.numeric(snakemake@params[['padj']])
gtf <- snakemake@params[['annotation']]
quant_folder <- snakemake@params[['quant_folder']]
comparison_condition_1 <- snakemake@params[['comparison_condition_1']]
comparison_condition_2 <- snakemake@params[['comparison_condition_2']]
conditions <- snakemake@params[['conditions']]

# List all directories in quant_folder
all_dirs <- list.dirs(quant_folder, recursive = FALSE, full.names = FALSE)
# Select directories based on the provided conditions and compare variables
# this requires the order of conditions match files but is more reliable than expecting all files to have names containing the condition name
samples_A <- file.path(quant_folder, all_dirs[which(grepl(comparison_condition_1, conditions))])
samples_B <- file.path(quant_folder, all_dirs[which(grepl(comparison_condition_2, conditions))])

# Make our salmon output into compatible kallisto type output
prepare_fish_for_sleuth(samples_A)
prepare_fish_for_sleuth(samples_B)

# Get transcript-gene id matches from gtf and ensure no target_id is empty
transcript_gene_idmatch <- gtf2ids(gtf)
transcript_gene_idmatch <- transcript_gene_idmatch %>% filter(target_id != "", 
                                                              target_id != "trnfM(CAT)", 
                                                              target_id != "trnS(TGA)", 
                                                              target_id != "trnY(GTA)")

# Calculate length- & library-normalised abundances. 
# Scale them to 1M reads for TPM values.
normalized_data <- fish4rodents(A_paths=samples_A, B_paths=samples_B, 
                         annot=transcript_gene_idmatch, scaleto=100000000)
print("Done normalizing DTU data")
# Call DTUs based on normalised data
called_DTUs <- call_DTU(annot=transcript_gene_idmatch, boot_data_A=normalized_data$boot_data_A, 
                    boot_data_B=normalized_data$boot_data_B, verbose=FALSE, abund_thresh=5, p_thresh=padj, dprop_thresh=cutoff, threads=2)
print("Done calling DTUs")
# Plot significance VS effect size:
#plot_overview(called_DTUs)

# Get all gene and transcript identifiers per category 
# (significant DTU, no DTU, Not Applicable):
DTU_ids <- get_dtu_ids(called_DTUs)
# Get all gene and transcript identifiers implicated in isoform switching:
DTU_genes <- DTU_ids[["DTU genes (gene test)"]]
writeLines(DTU_genes, snakemake@output[["DTU_ids_tsv"]])
print("Done writing DTU_ids.txt")

# Get all gene-wise results
write_tsv(as_tibble(called_DTUs$Genes), snakemake@output[["DTU_genes_tsv"]])
print("Done writing gene-wise DTU table")
# Get all transcript-wise results
write_tsv(as_tibble(called_DTUs$Transcripts), snakemake@output[["DTU_transcripts_tsv"]])
print("Done writing transcript-wise DTU table")

# Get abundance values
abundances1 <- called_DTUs$Abundances[[1]]
colnames(abundances1)[1:length(samples_A)] <- samples_A
abundances2 <- called_DTUs$Abundances[[2]]
colnames(abundances1)[1:length(samples_B)] <- samples_B
write_tsv(abundances1, snakemake@output[["abundance_cond1"]])
write_tsv(abundances2, snakemake@output[["abundance_cond2"]])
print("Done writing abundance tables")
