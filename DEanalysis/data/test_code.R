## Gene-level differential expression analysis using DESeq2

## Setup
### Bioconductor and CRAN libraries used
install.packages("tidyverse")


library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

## Load in data
data <- read.table("D:/Downloads/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("D:/Downloads/Mov10_full_meta.txt", header=T, row.names=1)

### Check classes of the data we just brought in
class(meta)
class(data)

View(meta)
View(data)

install.packages("ggplot2")
library(ggplot2)

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) + 
  xlim(-5, 500)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")

mean_counts <- apply(data[, 3:5], 1, mean)
variance_counts <- apply(data[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

View(counts(dds))

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="D:/Downloads/DEanalysis/data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA (for PC1 and PC2 only)
plotPCA(rld, intgroup="sampletype")

# Input is a matrix of log transformed values
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
pheatmap(rld_cor)

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run analysis
dds <- DESeq(dds)

## Plot dispersion estimates
plotDispEsts(dds)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")

library(apeglm)

## Control Vs Overexpression

## Define contrasts, extract results table, and shrink the log2 fold changes

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

resultsNames(dds)

res_tableOE <- lfcShrink(dds, coef=3, type="apeglm")

## MA Plot

plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))

class(res_tableOE)
mcols(res_tableOE, use.names=T)
library(dplyr)  
res_tableOE %>% data.frame() %>% View()

## Summarizing results
summary(res_tableOE)

## Extracting significant differentially expressed genes

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

library(tibble)
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sigOE

##Visualizing the results

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEGreport")

# load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport) 
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

# Create tibbles including row names
mov10_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Plot expression for single gene
plotCounts(dds, gene="MOV10", intgroup="sampletype") 

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="MOV10", intgroup="sampletype", returnData=TRUE)

# Plotting the MOV10 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))

## Order results by padj values
top20_sigOE_genes <- res_tableOE_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

## normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
  filter(gene %in% top20_sigOE_genes)

# Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:9], key = "samplename", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top20_sigOE)

gathered_top20_sigOE <- inner_join(mov10_meta, gathered_top20_sigOE)

## plot using ggplot2
ggplot(gathered_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

## Heatmap
### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names
norm_OEsig <- normalized_counts[,c(1,4:9)] %>% 
  filter(gene %in% sigOE$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### Annotate our heatmap (optional)
annotation <- mov10_meta %>% 
  select(samplename, sampletype) %>% 
  data.frame(row.names = "samplename")

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction

res_tableOE_tb <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

## Create a column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")

res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$gene[1:10]

View(res_tableOE_tb)

ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

sessionInfo()
