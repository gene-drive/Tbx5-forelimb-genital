########## Count matrix Heatmaps for RNA-seq data - last updated April 2023 by Aaron Alcala
# Run this in R/ RStudio

# Install packages: pheatmap, ComplexHeatmap, and InteractiveComplexHeatmap
# For more info on interactive heatmaps: https://github.com/jokergoo/InteractiveComplexHeatmap

install.packages("pheatmap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("InteractiveComplexHeatmap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEFormats")



# Load libraries
library("pheatmap")
library("ComplexHeatmap")
library("InteractiveComplexHeatmap")



##### Regular heatmap instructions

# Import your normalized count file from DEseq2. The column names are your replicates and the rows are each gene in the genome
dds_norm12a <- read.csv("DEseq2_normalized_counts_FL_2022_count_matrix12a_v2.csv", header=T, row.names=1)

# Import your Results Table file from DEseq2. The columns are baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
res12a_map <- read.csv("DEseq2_padj_0.05_FL_2022_count_matrix12a.csv")

# Check the stucture of these objects
head(dds_norm12a,5)
head(res12a_map,5)

# I renamed columns of res12a_map
colnames(res12a_map) <- c("row", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
head(res12a_map,1)

# I then renamed columns of dds_norm12a
colnames(dds_norm12a) <- c("ctrl_4", "ctrl_5", "ctrl_6", "ctrl_1", "ctrl_2", "ctrl_3", "KO_1", "KO_2", "KO_3", "KO_4", "KO_5", "KO_6")
head(dds_norm12a,1)

# I reordered the columns so my replicates were in proper order; now my control samples are the first six columns, and my mutants are the last six columns.
dds_norm12a <- dds_norm12a[, c(4,5,6,1,2,3,7,8,9,10,11,12)]

### Below I am sorting based on how I want the heatmap to appear. 
# Sort Results Table by padj
res12a_map_sorted <- res12a_map[order(res12a_map$padj),]
head(res12a_map_sorted, 4)

# Now I will filter so I only keep rows where padj value is less than 0.05
res12a_map_sorted2 <- subset(res12a_map_sorted, padj<0.05) 
head(res12a_map_sorted2, 4)

# Sort Results Table by log2FoldChange (ascending)
res12a_map_sorted3 <- res12a_map_sorted2[order(res12a_map_sorted2$log2FoldChange),]
head(res12a_map_sorted3, 4)


head(dds_norm12a, 4)

### Now I filter rows of the normalized count matrix to only show genes I'm interested in (based on the adjustments to the Results Table I did preceeding this section)

# In other words, the below command will only keep rows of the dds_norm12a count matrix that contain genes in res12a_map_sorted3
dds_norm12a_1 <- dds_norm12a[row.names(dds_norm12a) %in% res12a_map_sorted3$row, ]


class(dds_norm12a)
# [1] "data.frame"


# Sort dds_norm12a_1 so that gene IDs (rownames) are in same order as gene order set in res12a_map_sorted3 (The column name of the genes column is "row")
dds_norm12a_2 <- dds_norm12a_1[match(res12a_map_sorted3$row, rownames(dds_norm12a_1)), ]


# Save the newly filtered count matrix as as csv in case you want to use it later.
write.csv(as.data.frame(dds_norm12a_2), file="dds_norm12a_2.csv", row.names = TRUE)


# Convert data.frame into matrix
dds_norm12a_2 <- data.matrix(dds_norm12a_2, rownames.force = NA)


# THIS WILL MAKE THE BEAUTIFUL HEATMAP
dds_norm12a_2_map <- pheatmap(dds_norm12a_2, color=colorRampPalette(c("#4575B4", "white", "#E45139"))(25), cluster_rows=FALSE, show_rownames=FALSE, show_colnames=TRUE, cluster_cols=FALSE, scale="row")

dds_norm12a_2_map

###  Click the Export button above the heatmap image to save (This is in the Plots panel). 

# Saving as PDF gives good quality and is fastest way to export for basic presentations or lab meetings. Make sure to change size and orientation to your preference.

### To save for publication quality image, I save as a vector (NOT a raster image such as jpg, gif, etc.). 
# For some reason exporting as a .svg messes with the text. Therefore, I saved as .eps instead. Alternatively, you could save as a PDF (this is editable in Illustrator).
# To do this, click Export -> Save as Image. Then choose EPS as format. Then Resize the window so you get the desired image ratio for your heatmap. Then click save and it will save to your Working Directory.
# You can open the EPS file in Adobe Illustrator and adjust further. You can import the vector directly into Indesign or other programs like powerpoint. Alternatively, you can export your heatmap from Illustrator as a raster (jpg, tiff, or whatever you want).






##### Interactive heatmap instructions

# Make sure to read about InteractiveComplexHeatmap before starting: https://github.com/jokergoo/InteractiveComplexHeatmap
# Also: https://jokergoo.github.io/InteractiveComplexHeatmap/articles/InteractiveComplexHeatmap.html

# Load libraries if you haven't
library("ComplexHeatmap")
library("InteractiveComplexHeatmap")

# Viewing several examples by using below commands
htShinyExample(1.4)
htShinyExample(1.5)

# Because my heatmap was made using pheatmap, I need to use ComplexHeatmap:: to visualize it
dds_norm12a_2_map <- ComplexHeatmap::pheatmap(dds_norm12a_2, color=colorRampPalette(c("#4575B4", "white", "#E45139"))(25), cluster_rows=FALSE, show_rownames=FALSE, show_colnames=TRUE, cluster_cols=FALSE, scale="row")

dds_norm12a_2_map


# "Updating by draw() speeds up loading the Shiny application because draw() applies clusterings which is normally the most time-consuming step in heatmap generation."
dds_norm12a_2_map_draw <- draw(dds_norm12a_2_map)



# THIS WILL MAKE THE FANCY SCHMANCY INTERACTIVE HEATMAP
htShiny(dds_norm12a_2_map_draw)




# There are many more fancy ways you can split your data. Read the github for more info. 
