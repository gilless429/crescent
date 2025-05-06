# Load required packages
library(DESeq2, quietly = T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(ggrepel, quietly = T)
# Defining working directory (counts folder)
COUNTS_DIRECTORY <- snakemake@params[["input_path"]]
setwd(file.path(COUNTS_DIRECTORY))
# Define output path and create it if necessary
out_path <- snakemake@params[["output_path"]]
dir.create(out_path, showWarnings = FALSE, recursive = T)
# Pairs to compare together
PCA_labels <-  if (snakemake@params[["PCA_labels"]] == "yes") T else F
compare <- snakemake@params[["compare"]]
not_control <- compare[1]
control <- compare[2]

colData <- read.delim(snakemake@input[["colData"]], row.names=1, sep = '\t')
colData$condition <- as.factor(colData$condition)

if (snakemake@params[["quantif"]] == "salmon") {
    txi <- readRDS(file.path(COUNTS_DIRECTORY, "tximp.RDS"))
    dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
} else {
    # Grab the joined counts and colData matrixes made earlier by the count_matrix rule
    tsv_counts <- read_delim(snakemake@input[["count_matrix"]], delim ='\t')
    # Using 1st col (gene IDs) as row names as required for DESeq2, with transformation into dataframe as tibbles don't allow row names
    count_matrix <- as.data.frame(tsv_counts)[,2:ncol(tsv_counts)]
    row.names(count_matrix) <- tsv_counts[[1]]
    # Rounding and making into integers so if the input is fractional counts, it can still be read by DEseq2
    count_matrix <- round(count_matrix)
    count_matrix[,-1] <- lapply(count_matrix[,-1], as.integer)
    # Make sure our count_matrix has the same colnames as are used as rownames in colData
    colnames(count_matrix) <- rownames(colData)
    # Create DESeq2 object from our matrix and colData
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ condition)
    # Definition of reference point for each comparison
    dds$condition <- relevel(dds$condition, ref=control)
}
# Apply DESeq (normalization, DE analysis...)
dds <- DESeq(dds)
# Get results (lfc_shrink and not lfc_shrink)
res <- results(dds, contrast=c('condition', not_control, control), alpha=snakemake@params[["padj"]])
# resShr <- lfcShrink(dds, contrast=c('condition', not_control, control), res=res, type="ashr")
# Normalizing counts & outputting them
vsd <- varianceStabilizingTransformation(dds)
setwd(COUNTS_DIRECTORY)
write_delim(as_tibble(counts(dds, normalized=T), rownames = 'id'), paste0(out_path, '/normalized_counts_matrix.tsv'), delim = '\t')
# Creating and outputting a background of all the genes at least somewhat expressed in at least one sample
background <- rownames(filter_all(as.data.frame(counts(dds, normalized=T)), any_vars(. > 10)))
background <- unique(background)  # Making sure there is no duplicate
fwrite(list(background), file = paste0(out_path, "/background"))
# PCA plot of normalized counts
if (PCA_labels) {
    PCAplot <- plotPCA(vsd) + geom_label_repel(aes(label = colnames(vsd)), max.overlaps = Inf)
} else {
    PCAplot <- plotPCA(vsd)
}
ggsave(paste0(out_path, '/PCAplot.pdf'), plot=PCAplot)
ggsave(paste0(out_path, '/PCAplot.png'), plot=PCAplot)
# List of DEGs with user provided p-val and log2FC thresholds
resFilt <- res[ which(res$padj < snakemake@params[["padj"]] & abs(res$log2FoldChange) > snakemake@params[["log2FC"]]), ]
# Output tsv file with list of differentially expressed genes between the two concerned conditions
write_delim(as_tibble(resFilt, rownames="id"), paste(out_path, paste(not_control, 'VS', control, paste0("DEGstable_padj", snakemake@params[["padj"]],"_logFC", snakemake@params[["log2FC"]], ".tsv"), sep="_"), sep="/"), delim="\t")
