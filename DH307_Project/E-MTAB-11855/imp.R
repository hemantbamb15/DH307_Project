# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)  # For generating heatmaps

# Set working directory
setwd('C:/Users/heman/Desktop/E-MTAB-11855')

# Load count data
countData <- read.csv('gene_count_matrix.csv', header = TRUE, sep = ",")
head(countData)

# Load metadata
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
head(metaData)

# Construct DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ dex, tidy = TRUE)

# Run DESeq function
dds <- DESeq(dds)

# Check result table
res <- results(dds)
head(res)

# Set the output folder
output_folder <- "output_1"
dir.create(output_folder, showWarnings = FALSE)  # Create output folder

#Save the result table in csv format
write.csv(res, file.path(output_folder, "deseq_results_original.csv"), row.names = FALSE)

# Sort summary list with adjusted p-value
res <- res[order(res$padj),]
head(res)

#Save the result table new res(increasing p value) in csv format
write.csv(res, file.path(output_folder, "deseq_results_sorted.csv"), row.names = FALSE)



# Save all plots in a single PDF file
pdf_file <- file.path(output_folder, "all_plots.pdf")
pdf(pdf_file)

# Dispersion plot
plotDispEsts(dds)
bg_colors <- ifelse(metaData$dex == "untrt", "lightblue", "lightcoral")
# Original boxplots
par(mfrow=c(2, 2))
plotCounts(dds, gene = "ENSG00000272647.3|ENSG00000272647", intgroup = "dex",bg=bg_colors)
plotCounts(dds, gene = "ENSG00000187664.9|HAPLN4", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000260836.2|ENSG00000260836", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000249035.7|CLMAT3", intgroup = "dex")
par(mfrow=c(2, 2))
plotCounts(dds, gene = "ENSG00000254692.1|ENSG00000254692", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000287725.2|ENSG00000287725", intgroup = "dex")


# MA plot
plotMA(dds, ylim = c(-2, 2), main = "MA Plot")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Heatmap
counts_vsd <- assay(vsd)
row_means <- rowMeans(counts_vsd)
heatmap_data <- counts_vsd[row_means > 1, ]
heatmap_data <- heatmap_data[sample(nrow(heatmap_data), 100), ]

# Scale the data
heatmap_data_scaled <- scale(heatmap_data)

# Create the heatmap
pheatmap(heatmap_data_scaled, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE)

# Heatmap of log-transformed normalized counts. We will use the top 5 genes.
# Top 5 genes by p-value
top_hits <- res[order(res$padj), ][1:5, ]
top_hits <- row.names(top_hits)

# rlog transformation
rld <- rlog(dds, blind = FALSE)
pheatmap(assay(rld)[top_hits,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)
pheatmap(assay(rld)[top_hits,])


# PCA plot
plotPCA(vsd, intgroup = "dex")

dev.off()  # Close the PDF device

# Save results table
write.csv(res, file.path(output_folder, "deseq_results.csv"), row.names = FALSE)