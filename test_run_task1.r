#run test for a work at Dr. Wang Lab.
#task 1 on April 10
#Biosample (SAMN18442667 AND SAMN18442665) vs. (SAMN18442666 AND SAMN18442664) 
# to compare differentially-expressed genes and gene ontology enrichment analysis

#step 1: Retrieve RNA-seq Data (from NCBI)
# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SRAdb", "GEOquery"))
BiocManager::install("sra")

library(SRAdb)
library(RSQLite)

# Downloading SRA metadata SQLite file (It is 136 GB file)
sqlfile <- "C:/Users/SRAmetadb.sqlite/SRAmetadb.sqlite"
sra_con <- dbConnect(SQLite(), sqlfile)
dbListTables(sra_con)

# BioSamples clearly listed
# biosamples <- c("SAMN18442667", "SAMN18442665", "SAMN18442666", "SAMN18442664")
# Use SRA Experiment or Sample Accessions (NOT BioSample)
experiments <- c("SRX10431828", "SRX10431826", "SRX10431827", "SRX10431825")

conversion_results <- data.frame()

for (experiment in experiments) {
  result <- getSRA(search_terms = experiment, 
                   out_types = c("experiment", "run"),
                   sra_con = sra_con)
  
  conversion_results <- rbind(conversion_results, result)
}

# Final result with proper columns
final_results <- conversion_results[, c("experiment", "run")]
print(final_results)

#Step 1.2: Saving FASTQ files into a specific folder
# installing SRA Toolkit
# To download .sra files from NCBI
# To convert .sra to .fastq files
#Bash code:
# fasterq-dump --version
# prefetch SRR14055933
# fasterq-dump --split-files SRR14055933
# fasterq-dump --split-files --outdir "C:\Users\Tao\OneDrive\Desktop\ubc_dr_wang_lab\test" SRR14055933
# then, SRR14055935, SRR14055934, SRR14055936

#Step 2:Quality Control (QC) and Trimming
#Step 2.1: Run quality checks on your raw reads using FastQC
#Download FastQC for Windows.
#Add to PATH
#Java is needed. Thus, it should be installed before using.
#Run FastQC within one line from your updated folder "C:\Program Files\FastQC\".
#  cd "C:\Program Files\FastQC"
#  "C:\Program Files\FastQC\run_fastqc.bat" --outdir "E:\UBC_wang_qn1\fastqc_reports" "E:\UBC_wang_qn1\SRR14055933_1.fastq" "E:\UBC_wang_qn1\SRR14055933_2.fastq"
#  "C:\Program Files\FastQC\run_fastqc.bat" --outdir "E:\UBC_wang_qn1\fastqc_reports" "E:\UBC_wang_qn1\SRR14055934_1.fastq" "E:\UBC_wang_qn1\SRR14055934_2.fastq"
#  "C:\Program Files\FastQC\run_fastqc.bat" --outdir "E:\UBC_wang_qn1\fastqc_reports" "E:\UBC_wang_qn1\SRR14055935_1.fastq" "E:\UBC_wang_qn1\SRR14055935_2.fastq"
#  "C:\Program Files\FastQC\run_fastqc.bat" --outdir "E:\UBC_wang_qn1\fastqc_reports" "E:\UBC_wang_qn1\SRR14055936_1.fastq" "E:\UBC_wang_qn1\SRR14055936_2.fastq"
# thus, html and zip results are generated.

#Step 2.2: Trim adapters and low-quality bases using Trimmomatic
#Download Trimmomatic
#to check my installation
# dir "C:\my_program\Trimmomatic-0.39\Trimmomatic-0.39\adapters\TruSeq3-PE.fa"
#Trimmomatic Command in a single line
# cd /d C:\my_program\Trimmomatic-0.39\Trimmomatic-0.39
# since my command code with ^ does not work well, I write the command in a single line.
# java -jar trimmomatic-0.39.jar PE -threads 2 "E:\UBC_wang_qn1\SRR14055933_1.fastq" "E:\UBC_wang_qn1\SRR14055933_2.fastq" "E:\UBC_wang_qn1\SRR14055933_1_paired.fastq" "E:\UBC_wang_qn1\SRR14055933_1_unpaired.fastq" "E:\UBC_wang_qn1\SRR14055933_2_paired.fastq" "E:\UBC_wang_qn1\SRR14055933_2_unpaired.fastq" ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
# java -jar trimmomatic-0.39.jar PE -threads 2 "E:\UBC_wang_qn1\SRR14055934_1.fastq" "E:\UBC_wang_qn1\SRR14055934_2.fastq" "E:\UBC_wang_qn1\SRR14055934_1_paired.fastq" "E:\UBC_wang_qn1\SRR14055934_1_unpaired.fastq" "E:\UBC_wang_qn1\SRR14055934_2_paired.fastq" "E:\UBC_wang_qn1\SRR14055934_2_unpaired.fastq" ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
# similar for SRR14055935, SRR14055936

#Step 3: Mapping reads to the Reference Genome
#Step 3.1: Install HISAT2 from the website "https://www.di.fc.ul.pt/~afalcao/hisat2_windows.html"
#since there is no exe file in  this software, I have to add 3 files in the software by using Notepad: run_hisat2.bat, hisat2.bat, hisat2-build.bat.
# hisat2 -p 4 -x E:\ref_genomes\hisat2_mouse\mm10\genome -1 E:\UBC_wang_qn1\SRR14055933_1_paired.fastq -2 E:\UBC_wang_qn1\SRR14055933_2_paired.fastq -S E:\UBC_wang_qn1\SRR14055933.sam --summary-file E:\UBC_wang_qn1\SRR14055933_summary.txt

#Step 3.2: Convert, Sort, and Index the SAM File
#Install SAMtools
#However, no SAMtools for windows can be found.
#Ubuntu is installed.Then, by using Ubuntu, the code is below:
#sudo apt update
#sudo apt install samtools -y
#cd /mnt/e/UBC_wang_qn1
#samtools view -@ 4 -bS SRR14055933.sam > SRR14055933.bam
#samtools sort -@ 4 -o SRR14055933.sorted.bam SRR14055933.bam
#samtools index SRR14055933.sorted.bam
# I change to Ubuntu to get the SAM file as I did before with Windows CMD:
#/usr/bin/hisat2 -p 4 -x /mnt/e/ref_genomes/hisat2_mouse/mm10/genome -1 SRR14055935_1_paired.fastq -2 SRR14055935_2_paired.fastq -S SRR14055935.sam --summary-file SRR14055935_summary.txt
#Then I do bam, sort and index by using Ubuntu as above.

#Step 4: Gene-level quantification using featureCounts
#step 4.1: install featureCounts via Subread by using Ubuntu terminal.
# sudo apt install subread -y
#Step 4.2: Download Mouse Annotation File
# cd /mnt/e/ref_genomes/hisat2_mouse/mm10

# wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
# gunzip Mus_musculus.GRCm38.102.gtf.gz

#step 4.3: Run featureCounts
#cd /mnt/e/UBC_wang_qn1

#featureCounts -T 4 -a /mnt/e/ref_genomes/hisat2_mouse/mm10/Mus_musculus.GRCm38.102.gtf -o gene_counts.txt SRR14055933.sorted.bam SRR14055935.sorted.bam SRR14055934.sorted.bam SRR14055936.sorted.bam

#Step 5: Differential Expression Analysis in R (using DESeq2)
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")  # For shrinkage

# step 5.1:Read count matrix
library(DESeq2)

# Read the raw counts file (skip first comment lines)
counts <- read.table("E:/UBC_wang_qn1/gene_counts.txt", 
                     sep = "\t", header = TRUE, comment.char = "#", 
                     row.names = 1)

# Remove extra columns (e.g., "Chr", "Start", "End", etc.)
counts <- counts[,6:ncol(counts)]
colnames(counts) <- c("KO_1", "KO_2", "CTRL_1", "CTRL_2")

#step 5.2: Create DESeq2 dataset and run DESeq2
group <- factor(c("KO", "KO", "CTRL", "CTRL"))
coldata <- data.frame(row.names = colnames(counts), condition = group)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)
res <- lfcShrink(dds, coef = "condition_KO_vs_CTRL", type = "apeglm")  # Shrink fold changes
summary(res)

#5.3: Extract result

res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), "E:/UBC_wang_qn1/DEGs_KO_vs_CTRL.csv")

#Step 6: Visualization of Differential Expression
#install.packages("ggplot2")
#install.packages("pheatmap")

library(DESeq2)
library(ggplot2)
library(pheatmap)

#Step 6.1: Volcano Plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), size = 1.2) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

library(ggplot2)
library(ggrepel)

# Convert DESeq2 results to a dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Define significant DEGs (for labeling)
res_df$significant <- with(res_df, padj < 0.45 & abs(log2FoldChange) > 1)

# Select top genes to label: top 10 by smallest padj
top_genes <- res_df[res_df$significant, ]
top_genes <- top_genes[order(top_genes$padj), ][1:50, ]

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), size = 1.2) +
  scale_color_manual(values = c("gray", "red"), 
                     labels = c("Not Significant", "Significant (padj<0.45 & |LFC|>1)")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3.5, max.overlaps = 20) +
  labs(title = "Volcano Plot",
       subtitle = "Red = Significant DEGs (padj < 0.45 & |log2FC| > 1)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Gene Significance") +
  theme_minimal()

print(volcano_plot)

#Step 6.2: MA Plot
plotMA(res, ylim = c(-5, 5), main = "MA Plot")

#Step 6.3: PCA Plot (from normalized data)
library(DESeq2)
library(ggplot2)

# Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)

# Extract PCA data
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Make the PCA plot with tighter axis range
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  coord_cartesian(
    ylim = range(pcaData$PC2) + c(-1, 1),  # tighten Y axis range
    xlim = range(pcaData$PC1) + c(-1, 1)   # tighten X axis range
  ) +
  labs(
    title = "PCA of Samples (VST normalized)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(pca_plot)

# Step 6.4: Heatmap of Top DEGs
# Select top 30 DEGs
top_genes <- head(order(res$padj), 30)

# Extract normalized counts
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))  # Z-score scaling

# Annotate samples by group
annotation <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

# Plot heatmap
pheatmap(mat, annotation_col = annotation, show_rownames = TRUE, main = "Top 30 DEGs Heatmap")

#Step 7: GO Enrichment with clusterProfiler
#step 7.1: after reading the plots above, I want to save them for my report.
#volcano_plot
png("E:/UBC_wang_qn1/volcano_plot.png", width = 7, height = 6)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), size = 1.2) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

dev.off()

#PCA plot
png("E:/UBC_wang_qn1/pca_plot.png", width = 800, height = 600)

plotPCA(vsd, intgroup = "condition") + ggtitle("PCA of Samples (VST normalized)")

dev.off()

#MA plot
png("E:/UBC_wang_qn1/MA_plot.png", width = 800, height = 600)

plotMA(res, ylim = c(-5, 5), main = "MA Plot")

dev.off()

#Heatmap Plot
png("E:/UBC_wang_qn1/heatmap_top30.png", width = 900, height = 1200)

pheatmap(mat, annotation_col = annotation, show_rownames = TRUE, main = "Top 30 DEGs Heatmap")

dev.off()

#Step 7.2: Filter DEGs Based on Custom Thresholds
#Define DEGs (e.g., padj < 0.05 and |log2FoldChange| > 1)
# deg <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1), ]
deg <- res_df[which(res_df$padj < 0.4 & abs(res_df$log2FoldChange) > 0.1), ]
nrow(deg)
resLFC <- lfcShrink(dds, coef = "condition_KO_vs_CTRL", type = "apeglm")
summary(resLFC)
#deg_explore <- res_df[which(res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 1), ]
deg_explore <- res_df[which(res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 0.5), ]

nrow(deg_explore)

write.csv(deg_explore, "E:/UBC_wang_qn1/filtered_DEGs.csv")
###Due to limited statistical power (n=2 per group), no genes passed the FDR threshold. 
###However, exploratory analysis using unadjusted p-values (p < 0.05 & |log2FC| > 1) identified 8 genes. 
### exploratory analysis using unadjusted p-values (p < 0.05 & |log2FC| > 0.5) identified 10 genes. 
###These were used for preliminary GO enrichment.

###
##### I have to re-treat the GO enrichment because of deg data.
# Clean gene names
gene_symbols <- rownames(deg_explore)
gene_symbols <- sub("\\..*$", "", gene_symbols)

# Map to Entrez IDs
library(clusterProfiler)
library(org.Mm.eg.db)

gene_df <- bitr(gene_symbols,
                fromType = "ENSEMBL",
                toType   = "ENTREZID",
                OrgDb    = org.Mm.eg.db)

entrez_genes <- gene_df$ENTREZID
length(entrez_genes)  # Should return something like 6-8

ego <- enrichGO(gene         = entrez_genes,
                OrgDb        = org.Mm.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
head(ego)

print(barplot(ego, showCategory = 10, title = "Top GO Terms (Exploratory)"))
print(dotplot(ego, showCategory = 10, title = "GO Dot Plot (Exploratory)"))
png("GO_barplot.png")
barplot(ego, showCategory = 10, title = "Top GO Terms (Exploratory)")
dev.off()

png("GO_dotplot.png", width = 800, height = 600)
dotplot(ego, showCategory = 10, title = "GO Dot Plot (Exploratory)")
dev.off()

### My GO enrichment results are technically valid, but biologically very weak due to:
### A very small set of DEGs (only 7-10 genes used)
###Each GO term supported by just 1 gene
###Minimal statistical power (even though adjusted p-values were below 0.1)

###Use KEGG or Reactome Enrichment
library(clusterProfiler)
kegg <- enrichKEGG(gene = entrez_genes, organism = 'mmu')
dotplot(kegg)





# Step 7.3: GO Enrichment with clusterProfiler
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)
## Step 7.3.1: Prepare gene list
# Assume rownames(deg) are ENSEMBL IDs
gene_symbols <- rownames(deg)

# Convert to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_genes <- gene_df$ENTREZID

##Step 7.3.2: run go enrichment.
ego <- enrichGO(gene         = entrez_genes,
                OrgDb        = org.Mm.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",       # biological process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# Save results
write.csv(as.data.frame(ego), "E:/UBC_wang_qn1/GO_enrichment_BP.csv")

##step 7.3.3: Visualize GO terms (top 20)
summary(ego)

barplot(ego, showCategory = 20, title = "Top GO Terms (Biological Process)")
dotplot(ego, showCategory = 20, title = "GO Dot Plot")

## Use KEGG Enrichment
library(clusterProfiler)
library(org.Mm.eg.db)

# Subset DEG table with relaxed criteria
deg_explore <- res_df[which(res_df$pvalue < 0.2 & abs(res_df$log2FoldChange) > 0.3), ]

# Use ENSEMBL as the input type
gene_df <- bitr(deg_explore$gene,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# Extract unique ENTREZ IDs
entrez_genes <- unique(gene_df$ENTREZID)

# Check how many you got
length(entrez_genes)

##Step 8: Try gseKEGG() — a GSEA-style approach to rank all genes
##to try a different method to analyze further.
library(clusterProfiler)
library(org.Mm.eg.db)

# Step 8.1: Prepare gene list ranked by log2FoldChange
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$gene  # These are ENSEMBL IDs

# Step 8.2: Map ENSEMBL to ENTREZ
gene_map <- bitr(names(gene_list),
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

# Step 8.3: Merge gene_list and gene_map
gene_df <- data.frame(ENSEMBL = names(gene_list),
                      log2FC = gene_list)

merged_df <- merge(gene_df, gene_map, by = "ENSEMBL")

# Step 8.4: Remove duplicate ENTREZ IDs (keep the one with highest absolute LFC)
merged_df <- merged_df[!duplicated(merged_df$ENTREZID), ]

# Step 8.5: Prepare final named vector for GSEA
final_gene_list <- merged_df$log2FC
names(final_gene_list) <- merged_df$ENTREZID

# Step 8.6: Sort and remove NAs
final_gene_list <- sort(na.omit(final_gene_list), decreasing = TRUE)

# Step 8.7: Run GSEA (fgseaMultilevel will be automatically used)
gsea_kegg <- gseKEGG(geneList = final_gene_list,
                     organism = "mmu",
                     minGSSize = 10,
                     pvalueCutoff = 0.05)

# Step 8.8: Visualize (if results exist)
if (nrow(as.data.frame(gsea_kegg)) > 0) {
  dotplot(gsea_kegg, showCategory = 10, title = "GSEA KEGG Dotplot")
} else {
  print("No significant KEGG pathways enriched.")
}

write.csv(as.data.frame(gsea_kegg), file = "E:/UBC_wang_qn1/GSEA_KEGG_results.csv", row.names = FALSE)

#Step 9: generate the .tar.gz file (avoiding too large) by using Ubuntu:
# cd /mnt/e/UBC_wang_qn1
# tar --exclude='*.fastq' --exclude='*.sam' -czvf analysis_package.tar.gz *

