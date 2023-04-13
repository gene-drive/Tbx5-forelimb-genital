### RNA-seq analysis pipeline comparing Tbx5 cKO mutant forelimbs to controls - last updated April 2023
# Run this in R/ RStudio

# In this script, I initially used two separate control groups (heterozygotes and homozygous Tbx5 flox/flox-HoxB6Cre_negative) in the DESeq2 pipeline. Later in the script, I decide to remove the heterozygotes and only use the homozygous flox/flox as the control group to compare to Tbx5 cKO mutants.

# The annotation file used was: gencode.vM10.annotation.gtf



### For referennce and vignette, see: Analyzing RNA-seq data with DESeq2
# Michael I. Love, Simon Anders, and Wolfgang Huber
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

install.packages("Rsubread")
install.packages("DESeq2")

# Set working directory to where you want outputs to go
setwd("C:/Users/aaron/Dropbox/1-UGA/1-Research/Menke Lab/RNA-seq/FL_RNA-seq_2022")

# Below command lists files from directory where my alignment (.bam) files are for this analysis. I will use this object in the featureCounts command:
fls10a <- list.files( "E:/1-UGA-4 TB HDD/FL_RNAseq bams/renamed for DESeq2", pattern="bam$", full.names = TRUE)

fls10a

library("Rsubread")

# featureCounts command below!
print(Sys.time())
fc10a <- featureCounts(files=fls10a,
                       annot.ext="E:/1-UGA-4 TB HDD/gencode.vM10.annotation.gtf", 
                       isGTFAnnotationFile=TRUE,
                       GTF.featureType="exon",
                       GTF.attrType="gene_name",
                       allowMultiOverlap=FALSE,
                       isPairedEnd=FALSE,
                       strandSpecific=2)

# Above command takes about 1 minute per replicate on my desktop PC. Can take a coffee break.

# Extract count matrix from fc object
count_matrix10a <- as.matrix(fc10a$count)
head(count_matrix10a,2)
View(count_matrix10a)

########################

# Rename column names of the matrix:
colnames(count_matrix10a) <- c("ctrl_1", "ctrl_2", "ctrl_3", "ctrl_4", "ctrl_5", "ctrl_6", "ctrlb_10", "ctrlb_11", "ctrlb_12", "ctrlb_7", "ctrlb_8", "ctrlb_9", "KO_1", "KO_2", "KO_3", "KO_4", "KO_5", "KO_6")

head(count_matrix10a,2)
View(count_matrix10a)

#To save matrix as csv
write.csv(count_matrix10a, "FL-RNAseq_2022_count_matrix10a.csv")

#### Before proceeding, make a colData .csv file containing sample info. My colData_FL_allreps_v4.csv file has three columns that look like this:
# SampleName	Genotype	Replicate
# ctrl_1	ctrl	1
# ctrl_2	ctrl	2
# ctrl_3	ctrl	3
# ctrl_4	ctrl	4
# ctrl_5	ctrl	5
# ctrl_6	ctrl	6
# ctrlb_10	ctrlb	10
# ctrlb_11	ctrlb	11
# ctrlb_12	ctrlb	12
# ctrlb_7	ctrlb	7
# ctrlb_8	ctrlb	8
# ctrlb_9	ctrlb	9
# KO_1	KO	1
# KO_2	KO	2
# KO_3	KO	3
# KO_4	KO	4
# KO_5	KO	5
# KO_6	KO	6

#colData is your sampleinfo
# Load the csv into Rstudio
csvfile10a <- file.path("C:/Users/aaron/Dropbox/1-UGA/1-Research/Menke Lab/RNA-seq_scripts","colData_FL_allreps_v4.csv")

coldata10a <- read.csv("C:/Users/aaron/Dropbox/1-UGA/1-Research/Menke Lab/RNA-seq_scripts/colData_FL_allreps_v4.csv",row.names=1)
head(coldata10a,19)

# Check if naming is consistent between colData and count matrix
all(rownames(coldata10a) %in% colnames(count_matrix10a))
# [1] TRUE

all(rownames(coldata10a) == colnames(count_matrix10a))
# [1] TRUE

# Proceed to DGE analysis below!!!!!

library("DESeq2")
dds10a <- DESeqDataSetFromMatrix(countData = count_matrix10a, colData = coldata10a, design = ~ Genotype)

#I got this message:
# Warning message:
#   In DESeqDataSet(se, design = design, ignoreRank) :
#   some variables in design formula are characters, converting to factors

# Generate normalized count matrix
dds10a <- estimateSizeFactors(dds10a)
dds_norm10a <- counts(dds10a, normalized=TRUE)
write.csv(as.data.frame(dds_norm10a), file="DEseq2_normalized_counts_FL_2022_count_matrix10a.csv", row.names = TRUE)

# Make the reference (e.g. control, or untreated samples) as the first level for selected column name (e.g. genotype/tissue) from sampleinfo	
dds10a$Genotype <- relevel(dds10a$Genotype, "ctrl")

# Filter out genes with low numbers of reads, which have no information about the amount of gene expression
dds10a <- dds10a[rowSums(counts(dds10a)) > 1, ]

# Run the differential expression pipeline 
deseq_dds10a <- DESeq(dds10a)

# Build results table at selected significance level and sort them by selected column of the results table (e.g. padj)
res10a <- results(deseq_dds10a, alpha=0.05)
res_sorted10a <- res10a[order(res10a$padj),]

# Display the summary of the results	
summary(res_sorted10a)

# out of 29157 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 179, 0.61%
# LFC < 0 (down)     : 194, 0.67%
# outliers [1]       : 114, 0.39%
# low counts [2]     : 6172, 21%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save the results table for downstream analyses!	
write.csv(as.data.frame(res_sorted10a), file="DEseq2_padj_0.05_FL_2022_count_matrix10a.csv", row.names = TRUE)

# Done! 



# You can save your Workspace via the Environment tab (or using the below command). Then you can load this workspace up and pick up where you left off from last time with all your objects!
save.image("E:/1-UGA-4 TB HDD/FL TBX5 RNAseq analysis_workspace_Feb2022.RData")





###### To re-analyze my data with a subset of the replicates, I used the following commands to manipulate my matrix object. 

# Removing several replicates because I no longer want them used as control samples:
count_matrix12a <- count_matrix10a[,-c(1:6)]
count_matrix12a

# Removing a replicate because we determined it was an outlier:
count_matrix14a <- count_matrix12a[,-c(12)]
count_matrix14a

# After I trim out these replicates (columns) from my count_matrix, I used the matrix for the DESeq2 analysis.


# You can reorder matrix columns using:
dds_norm14a <- dds_norm14a[, c(4,5,6,1,2,3,7,8,9,10,11)]





########## QC and Visualization ##############

# Plotting PCA and Heatmap of the sample-to-sample distances. Good for finding outliers.

# For publication-quality images, export figures as PDF so you can easily manipulate colors, text, etc. in Adobe Illustrator. Usually, you do NOT want to save as a raster image such as jpg, gif, etc.
# Alternatively save your figure as .eps vector file. For some reason exporting as a .svg messes with the text.
# You can import the vector file directly into Adobe InDesign or other programs like PowerPoint.


install.packages("ggrepel")
install.packages("pheatmap")
library("pheatmap")
library("ggrepel")
library("ggplot2")

# This code is in the DESeq2 vignette
transformed_count_matrix10a <- vst(dds10a, blind=FALSE)
class(transformed_count_matrix10a)

# Plot the PCA!
plotPCA(transformed_count_matrix10a, intgroup=c("Genotype")) + geom_text_repel(aes(label=factor(colnames(transformed_count_matrix10a))), size = 3)


## Heatmap of the sample-to-sample distances
# A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. 

# Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists_10a <- dist(t(assay(transformed_count_matrix10a)))

# Below commands will generate the sample-to-sample heatmap
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists_10a)
rownames(sampleDistMatrix) <- paste(colnames(transformed_count_matrix10a))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists_10a,
         clustering_distance_cols=sampleDists_10a,
         col = colors)



####### Count matrix Heatmaps, Volcano plots, and Gene Ontology Enrichment Analyses ###########

# To visualize your (normalized) counts data as an interactive heatmap online, use Morpheus. Good for quickly looking at counts of genes of interest. 
# https://software.broadinstitute.org/morpheus/


# To generate a Heatmap of the count matrix, I used the pheatmap and ComplexHeatmap packages
# https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
# Please see this separate script I wrote for how I generated my Heatmaps:
# https://github.com/gene-drive/Tbx5-forelimb-genital/blob/307cb94428ff260fcce215ee8ce83b98c24b939d/Heatmaps_RNA-seq


# Customizable and easy-to-use webtool for Volcano Plots: ggVolcanoR
# https://ggvolcanor.erc.monash.edu/


# Instead of DAVID, I used the ShinyGO webtool for GO analyses and to generate pretty lollipop plots and KEGG network diagrams. Sometimes it loads slowly though.
# http://bioinformatics.sdstate.edu/go/


# Don't be afraid to email the authors/maintainers of the software as they are usually pretty responsive and even willing to squash out bugs or make the tools more user-friendly.
