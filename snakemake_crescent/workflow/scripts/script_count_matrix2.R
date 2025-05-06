# Import relevant libraries
suppressMessages(library(tidyverse))

# Set up empty matrix for later use and define the path in which count files are to be found
count_matrix <- matrix()
counts_path <- snakemake@params[["path"]]
# Define colData for use in DESeq2
sample_names <- snakemake@params[["sample_names"]]
colData <- tibble("samples" = sample_names, "condition" = snakemake@params[["conditions"]])
write_delim(colData, file.path(counts_path, "colData.tsv"), delim="\t")

if (snakemake@params[["quantif"]] == "salmon") {
  library(tximport)
  files <- file.path(list.dirs(path = counts_path, recursive = F), "quant.sf")
  names(files) <- sample_names
  # Importing tx2gene correspondance matrix (matching transcript ids to gene ids)
  tx2gene <- read_delim(file.path(snakemake@params[["working_dir"]], 'RNAseq', 'tx2gene.tsv'), delim = '\t')
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  count_matrix_sorted <- as_tibble(txi$counts, rownames = 'id')
  rownames(colData) <- colnames(txi$counts)
  # Saving the tximport object to RDS to use it 
  saveRDS(txi, file = file.path(counts_path, 'tximp.RDS'))
} else {
  # Define comment character and column containing counts depending on the quantification tool used
    if (snakemake@params[["quantif"]] == "featurecounts") {
      comments <- '#'
      count_column <- 7
    } else if (snakemake@params[["quantif"]] == "htseq") {
      comments = '__'
      count_column <- 2
    }
  # Read the count files and add them to the count matrix
  for (file in list.files(counts_path, pattern=".txt$", full.names=TRUE)) {
    table <- read_delim(file, delim = "\t", escape_double = FALSE, comment = comments, trim_ws = TRUE)
    count_matrix <- bind_cols(count_matrix, table[,count_column])
  }
  # Remove the first column of the count matrix, which is empty (due to how it was created)
  count_matrix <- count_matrix[,-1]
  
  # Define column names for the count matrix based on the file names for the counts (stripping away everything else to keep only the sample name)
  file_names <- list.files(counts_path, pattern=".txt$")
  if (length(sample_names) != ncol(count_matrix)) {
    sample_names <- unlist(lapply(file_names, function(x) strsplit(strsplit(x[[1]], '/')[[1]][[length(strsplit(x[[1]], '/')[[1]])]], snakemake@params[["sep_file"]])[[1]][[snakemake@params[["index_sample"]]]]))
  }
  colnames(count_matrix) <- sample_names
  
  # Define 'id' column, providing gene identification per row, then sort based on it to get rows in order of gene name and put 'id' column first
  count_matrix['id'] <- table[,1]
  count_matrix_sorted <- count_matrix[order(count_matrix[['id']]),]
  count_matrix_sorted <- count_matrix_sorted %>% select('id', everything())
}

# Making sure the id column comes first
count_matrix_sorted <- count_matrix_sorted %>% select('id', everything())
# Export the count matrix to a file, which will be used by the global report (& the differential expression analysis by DESeq2 unless salmon was used)
write_delim(count_matrix_sorted, file.path(counts_path, "count_matrix.tsv"), delim="\t")