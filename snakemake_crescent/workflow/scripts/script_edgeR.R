# Load required packages
library(edgeR, quietly = T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(ggrepel, quietly = T)
library(RColorBrewer)
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
# Definition of reference point for each comparison
colData$condition <- relevel(colData$condition, ref=control)

if (snakemake@params[["quantif"]] == "salmon") {
    txi <- readRDS(file.path(COUNTS_DIRECTORY, "tximp.RDS"))
    # Following standard edgeR tximport procedure
    cts <- txi$counts
    normMat <- txi$length
    normMat <- normMat/exp(rowMeans(log(normMat)))
    normCts <- cts/normMat
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)
    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)
    DGE_object <- DGEList(cts)
    DGE_object <- scaleOffset(DGE_object, normMat)
    # filtering using the design information
    design <- model.matrix(~condition, data = colData)
    keep <- filterByExpr(DGE_object, design)
    DGE_object_ftr <- DGE_object[keep, ]
} else {
    # Grab the joined counts and colData matrixes made earlier by the count_matrix rule
    tsv_counts <- read_delim(snakemake@input[["count_matrix"]], delim ='\t')
    # Using 1st col (gene IDs) as row names as required for DESeq2, with transformation into dataframe as tibbles don't allow row names
    count_matrix <- as.data.frame(tsv_counts)[,2:ncol(tsv_counts)]
    row.names(count_matrix) <- tsv_counts[[1]]
    # Rounding and making into integers so if the input is fractional counts, it can still be read by edgeR
    count_matrix <- round(count_matrix)
    count_matrix[,-1] <- lapply(count_matrix[,-1], as.integer)
    # Make sure our count_matrix has the same colnames as are used as rownames in colData
    colnames(count_matrix) <- rownames(colData)
    # Create initial edgeR object from our matrix and colData
    DGE_object<- DGEList(counts=count_matrix,group=colData$condition)
   # Low counts filtering
    design <- model.matrix(~condition, data = colData)
    keep <- filterByExpr(DGE_object, design)
    DGE_object_ftr <- DGE_object[keep, , keep.lib.sizes=FALSE]
}
# TMM normalisation
DGE_object_TMM <- calcNormFactors(DGE_object_ftr, method="TMM")
# Global dispersion estimate
DGE_object_DEA <- estimateDisp(DGE_object_TMM, design)
# Model fitting to count data
fit <- glmQLFit(DGE_object_DEA, design)
# Genewise hypothesis tests
coef_index <- which(colnames(design) == paste0("condition", not_control))
results <- glmQLFTest(fit, coef = coef_index)
# Overall results recovering
dea_res <- topTags(results, n = Inf, adjust.method = "BH")
# Output read counts (normalized TMM read counts)
write_delim(as_tibble(cpm(DGE_object_TMM), rownames = 'id'), paste0(out_path, '/normalized_counts_matrix.tsv'), delim = '\t')
# Counts transformation for PCA
log_cpm <- cpm(DGE_object_TMM, log=TRUE, prior.count=2)
# Creating and outputting a background of all the genes at least somewhat expressed in at least one sample
background <- names(keep[keep==TRUE])
background <- unique(background)  # Making sure there is no duplicate
fwrite(list(background), file = paste0(out_path, "/background"))
# PCA plot of transformed counts
pdf(paste0(out_path, "/PCAplot.pdf"))
plotMDS(log_cpm, gene.selection = "common", col=brewer.pal(min(9, length(levels(colData$condition))), "Set2")[colData$condition], plot=T)
dev.off()
png(paste0(out_path, "/PCAplot.png"), width = 700, height = 700)
plotMDS(log_cpm, gene.selection = "common",  col=brewer.pal(min(9, length(levels(colData$condition))), "Set2")[colData$condition], plot=T)
dev.off()

# List of DEGs with user provided p-val and log2FC thresholds
resFilt <- dea_res$table[ dea_res$table$FDR < snakemake@params[["padj"]] & abs(dea_res$table$logFC) > snakemake@params[["log2FC"]], ]
# Output tsv file with list of differentially expressed genes between the two concerned conditions
write_delim(as_tibble(resFilt, rownames="id"), paste(out_path, paste(not_control, 'VS', control, paste0("DEGstable_padj", snakemake@params[["padj"]],"_logFC", snakemake@params[["log2FC"]], ".tsv"), sep="_"), sep="/"), delim="\t")
