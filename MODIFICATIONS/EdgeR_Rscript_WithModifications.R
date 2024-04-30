#### Install the EdgeR package if you have not already
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

## Load the edgeR library 
library(limma)
library(edgeR)
library(dplyr)
library(stringr)
## Use the Session menu to set working directory To Source File Directory

# Reading in the feature count file as "counts.df"
countdata <- read.csv("gene_count_matrix.csv", row.names="gene_id")
# Clean the gene IDs (currently your row names)
clean_gene_ids <- str_replace(rownames(countdata), "gene-", "")
clean_gene_ids <- str_remove(clean_gene_ids, "\\|.+") 

# Update the row names with the cleaned versions
rownames(countdata) <- clean_gene_ids


# Printing the start of the counts.df object in R...
dim(countdata)
head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table(file.choose(), header=TRUE,row.names = 1))
dim(coldata)
head(coldata)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

# Subsetting gene counts according to experimental condition
counts_control.df  <- countdata[,c("SRR6825930", "SRR6825925", "SRR6825926", "SRR6825938", "SRR6825934", "SRR6825935")]
counts_treated.df <- countdata[,c("SRR6825928", "SRR6825932", "SRR6825933", "SRR6825940", "SRR6825941", "SRR6825939")]

# Printing the structure of the gene counts set and subsets
str(countdata)
str(counts_control.df)
str(counts_treated.df)


#If you dont want to use RSD test, then creating a DGEList object using the original unfiltered gene counts
counts.DGEList <- DGEList(counts = countdata, genes = rownames(countdata))

# Printing the design table
print(coldata)

# for unfiltered(with noRSD filtering) data use this
summary(colnames(countdata) == rownames(coldata))

# Add grouping information to DGEList object
counts.DGEList$samples$group <- as.factor(coldata$Treatment)

# Printing counts.DGEList
counts.DGEList
dim(counts.DGEList)

# Creating an object to filter genes with low expression
counts.keep <- filterByExpr(counts.DGEList, min.count=20)
#We want to change the min.count just like we changed the prefiltering parameter in Deseq2 pipeline, so this value can be adjusted to 0, 10, 50 and 100.
#I am keeping it to 20 here just because it is what we use in our standard Deseq2 pipeline)
#counts.keep <- filterByExpr(counts.DGEList, min.count=0)
#counts.keep <- filterByExpr(counts.DGEList, min.count=10)
#counts.keep <- filterByExpr(counts.DGEList, min.count=50)
#counts.keep <- filterByExpr(counts.DGEList, min.count=100)
summary(counts.keep)

# Filtering lowly expressed genes
filtered.DGEList <- counts.DGEList[counts.keep, , keep.lib.sizes = FALSE]
dim(filtered.DGEList)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in counts.keep
length(counts.keep[counts.keep == TRUE]) == dim(filtered.DGEList)[1]

# Printing the normalisation factors for the libraries
filtered.DGEList$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
norm.DGEList <- calcNormFactors(filtered.DGEList)
norm.DGEList$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition <- coldata$Treatment
disp.DGEList <- estimateDisp(norm.DGEList, design = model.matrix(~condition))
disp.DGEList

# Exact tests for differences between experimental conditions
std_treated.DGEExact <- exactTest(disp.DGEList, pair = c("Control",
                                                         "Heat_treated"))

# Extracting most differentially expressed genes from exact tests
std_treated.topTags <- topTags(std_treated.DGEExact,n = nrow(std_treated.DGEExact$table))

# Printing the most differentially expressed genes
std_treated.topTags

# extract significant differentially expressed genes, sort, & write to csv
resOrdered <-std_treated.topTags$table[order(std_treated.topTags$table$FDR),]
resSig <- resOrdered[resOrdered$FDR<0.05,]
print(resOrdered)

out <- resOrdered %>% dplyr::select(logFC, logCPM, PValue, FDR) %>% dplyr::rename(meanExpr = logCPM, pval = PValue, adj.pval = FDR)
write.csv(as.data.frame(out),file = paste0("C:/Users/arpan/OneDrive/Desktop/Final project/edgeR3.csv")) #write results to a new csv

#Calculating the total number of differentially expressed genes at FDR< 0.1
#Again, I am using 0.1 as the FDR value as this is what the standard pipeline in Deseq2 uses. We can change the p.value to 0.05 and 0.01 here but keep in mind that
#that while changing p.value we are keeping min.count in filterbyexpr step to 20 
de1 <- decideTestsDGE(std_treated.DGEExact , adjust.method="BH", p.value=0.1)
#de1 <- decideTestsDGE(std_treated.DGEExact , adjust.method="BH", p.value=0.05)
#de1 <- decideTestsDGE(std_treated.DGEExact , adjust.method="BH", p.value=0.01)
summary(de1)

# Create the MA-plot 
#Whatever p.value we use in above code, we use the same here
edgeR::plotMD.DGEExact(std_treated.DGEExact,p.value=0.1)
#edgeR::plotMD.DGEExact(std_treated.DGEExact,p.value=0.05)
#edgeR::plotMD.DGEExact(std_treated.DGEExact,p.value=0.01)

#..............heatmap construction.........................................
library("RColorBrewer")
library("pheatmap")

logcpm <- cpm(norm.DGEList, log=TRUE) 
## Making heatmap of count matrix for significantly expressed genes
#Prepare Expression Matrix by subsetting the log-CPM matrix to include only significant genes
topVarGenes <- order(apply(logcpm, 1, var), decreasing=TRUE)[1:50] 
expression_matrix <- logcpm[topVarGenes, ]

#expression_matrix <- logcpm[rownames(resSig), ]

# Choose Color Palette
color_palette <- brewer.pal(9, "RdYlBu")  # Adjust the palette name if desired

treatment_groups <-  coldata$Treatment 
annotation_df <- data.frame(Treatment = treatment_groups) 
rownames(annotation_df) <- colnames(expression_matrix)

# Generate Heatmap
pheatmap(expression_matrix,
         color = color_palette,
         scale = "row",          # Scale expression values by row
         cluster_rows = TRUE,    
         cluster_cols = TRUE,   
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap of Significant Genes",
         annotation_col = annotation_df) # Add a title
## If you want to make heatmap for all the genes after filtering steps, use logcpm instead of expression matrix
#pheatmap(expression_matrix,
#         color = color_palette,
#       cluster_rows = TRUE,    
#         cluster_cols = TRUE,  
#         show_rownames = TRUE, 
#        show_colnames = TRUE,
#         main = "Heatmap of Significant Genes",
#         annotation_col = annotation_df) # Add a title

#Heatmap for sample to sample distances

expression_matrix2 <- t(logcpm[rownames(resSig), ])
# Calculate sample distances
sampleDists <- dist(expression_matrix2, method = "euclidean") 
sampleDistMatrix <- as.matrix(sampleDists)

# Color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Heatmap of Sample-to-Sample Distances") 

df <- data.frame(coldata$Treatment, logcpm[topVarGenes,])

# Perform PCA using limma's plotMDS function
plotMDS(expression_matrix, main = "PCA Plot", col = as.numeric(coldata$Treatment)) 

#GSEA rank file making

DGE_Anno_Rank <-  within(std_treated.topTags$table, rank <- sign(logFC) * -log10(PValue))
DGE_Anno_Rank 

#subset the results so only Gene Name and rank
DGErank = subset(DGE_Anno_Rank, select = c(genes,rank) )
DGErank

#sebset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)



### IF NEEDED....remove the "gene-" from row names
## https://stackoverflow.com/questions/39897155/replace-and-remove-part-of-string-in-rownames/39897315  "URS000075AF9C-snoRNA_GTATGTGTGGACAGCACTGAGACTGAGTCT"    to   "snoRNA"
## We can use gsub to match one of more characters that are not a - ([^-]+) from the start (^) of the string followed by 
## a - or (|) one or more characters that are not an underscore ([^_]+) until the end of the string ($)  and replace it with blanks ("").
## gsub("^[^-]+-|_[^_]+$", "", rownames(df))

#gene <-gsub("^[^-]+-", "", rownames(DGErank))
#DGErankIDs  <-cbind(gene,DGErank)
#head(DGErankIDs)
#summary(DGErankIDs)  

#write.csv(as.data.frame(DGErankIDs), file="DGErankIDs.csv", row.names=FALSE)
#write.table(as.data.frame(DErank_withName),file="DEGrank_withName.txt",quote=FALSE, row.names=FALSE,sep="\t")
write.table(as.data.frame(DGErank_withName),file="DGErank_withName.rnk",quote=FALSE, row.names=FALSE,sep="\t")
##................................Normalization...................................................
# Create the data frame
NormTransExp <- data.frame(DGErank$genes, logcpm) 

# Export as CSV
write.csv(NormTransExp, file="NormTransExpressionData.txt", row.names=FALSE)