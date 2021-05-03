#### Load and install required packages ####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics") #For PLS-DA
BiocManager::install("tximport") #To import gene counts from Salmon
BiocManager::install("GenomeInfoDb")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
install.packages("geomorph") #For shape landmark alignment and 2BPLS
install.packages("pls")

library(mixOmics)
library(tximport)
library(GenomeInfoDb)
library(org.Dm.eg.db) 
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(geomorph)
library(pls)

#### Set up data ####

#set up the working directory
setwd("/Users/amandan/Desktop/grad_classes/bio722/individual_project/data")

### Load in gene counts 
#define path to count files
quant_files <- file.path("dgrp_counts", list.files("dgrp_counts"), "quant.sf")

#read in the sample metadata
sample_info <- read.csv("dgrp-sample_info.csv")

#order sample_info so that it matches quant_files
sample_info <- sample_info[order(sample_info$fastqFile),]

#add a column for "Line" so that lines can be matched with the shape data Line column
sample_info$LineNum <- gsub("_|_1|_2|-1|-2|_T1|T2", "", x = sample_info$lineName)


### Load in shape data
raw_shape <- read.csv("BothLabs_Wings_28Oct.csv")

#I am only going to work with the flies from the Dworkin lab, not the Houle lab, so subset these
shape_data <- raw_shape[grep("Dwo", raw_shape$Lab),]

#Combine $Line and $Ind into one column for easier identification
shape_data$ID <- paste(shape_data$Line, shape_data$Ind, sep = "_")

#Adjust csize (centroid size, a measure of wing size), needs to be multiplied by the scale
shape_data$Csize <- shape_data$Csize * shape_data$Scale

#There are 4 replicates for each line. Prior work has shown the differences between replicates
#(and sexes) are negligible, so to make things easy I am going to only use the 3rd replicate for each sample
shape_data <- shape_data[grep("rep3", shape_data$Rep),]

#### Data set matching ####
#Since the exact lines measured for the RNA seq and shape analysis might not be the same,
#I want to subset both datasets so that only individuals from matching lines are kept (otherwise, they are not really comparable)
as.integer(sample_info$LineNum)
shape_data$LineNum

gene_sample <- sample_info[sample_info$LineNum %in% shape_data$LineNum,]
dim(gene_sample) #only contains the 81 matching lines

shape_sample <- shape_data[shape_data$LineNum %in% sample_info$LineNum,]
dim(shape_sample) #now contains 4496 shape measurements

#Back to the counts, use the metadata to give names to the quant files, while removing file extensions from sample names 
names(quant_files) <- gsub(".txt.gz|.fastq.gz", "", x = sample_info$fastqFile) 
#subset quant_files so that it contains only the matching lines, as in gene_sample
quant_files <- quant_files[sample_info$LineNum %in% shape_data$LineNum]

#Now that the data has been subsetted and R knows which files to read in, I am going to import the sample counts
#Setting up my gene to transcript identifier 
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(x = txdb, keys = k, "GENEID", "TXNAME")

#Import the count files using Tximport
txi <- tximport(quant_files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance="scaledTPM")

cts <- txi$counts #Make the counts their own object
dim(cts) #Looks to be the right number of samples (87 columns)

#For a quick check to make sure things imported okay, I am going to check out the counts for 
#vestigial (vg), which is known to be expressed in the wing
cts[rownames(cts) == 'FBgn0003975',] #Looks good 

#checking eyeless FBgn0005558, should not be in the wing 
cts[rownames(cts) == 'FBgn0005558',] #Looks good (not exactly 0, which is to be expected, but close)

#To make my life easier, I am only going to look at a subset of genes, one for which we expect to be expressed during
#wing development in the wing tissue
#I downloaded a list of genes associated with wing development from FlyBase.org
wing_genes <- read.csv("FlyBase_IDs_wing_development.txt", header = FALSE)
cts_wing <- cts[rownames(cts) %in% wing_genes$V1,]
dim(cts_wing) #This leaves us with 372 genes from an original value of 13701

#Turn cts into a 3D array so that it behaves well with geomorph two.b.pls
gene_arr <- array(data = cts_wing, dim = c(372, 1, 87), dimnames = list(rownames(cts_wing), "count", gene_sample$lineName))

#The shape coordinates have already been aligned in order to remove the effects of variation
#in location, orientation, and scale. This ensures what we are left with is shape data
#and not size data
coords <- shape_sample[,10:105] #extract just the landmark coordinates 
coords_arr <- arrayspecs(coords, p = 48, k = 2) #turn the landmark coordinates into a 3D array, which behaves better with geomorph
plotAllSpecimens(coords_arr) #plot the landmark coordinates. Does this look like a wing to you?

#To finish setting up the data, I am going to rename the dimensions of coords_arr so things are easier to look at and keep track of
lm_num <- 1:48
dim_names <- c("X", "Y")
line_names <- as.character(shape_sample$ID)
dimnames(coords_arr) <- list(lm_num, dim_names, line_names)

#### 2BPLS Analysis ####
#geomorph two.b.pls
#The considerations here are that we need a 3d landmark array and a 2d variable matrix (for gene counts)
#Additionally, the function assumes that all specimens are in the same order
#Again, I am going to ignore replicates and just take one of each sample for both the 
#shape data and the expression data

#Turn the gene data into a 2d matrix (rows = specimens, columns = variables)
gene_mat <- two.d.array(gene_arr)
rownames(gene_mat)

gene_subset <- gene_mat[-grep("_2|-2|_T2", rownames(gene_mat)),] #remove 2nd reps
rownames(gene_subset) <- gsub("_1|-1|_T1", "", rownames(gene_subset)) #fix row names 
#Remove two (out of four) of the 440 samples manually
grep("440", rownames(gene_subset))
rownames(gene_subset)[grep("440", rownames(gene_subset))]
gene_subset <- gene_subset[-c(35,56,73),]
grep("440", rownames(gene_subset))

dim(gene_subset) #left with 72 unique samples


#Do the same for shape_sample, this time using match()
shape_subset <- shape_sample[match(rownames(gene_subset), shape_sample$LineNum),]
#Remake the array
shape_subset_arr <- arrayspecs(shape_subset[,10:105], p = 48, k = 2)
line_names <- as.character(shape_subset$LineNum)
dimnames(shape_subset_arr) <- list(lm_num, dim_names, line_names)


#two.b.pls, consider adding log size as another matrix
gm_2bpls <- two.b.pls(shape_subset_arr, gene_subset)
gm_2bpls$right.pls.vectors #first column is vecotr with greatest amount of covariation
gm_2bpls$left.pls.vectors #shape effects
str(gm_2bpls$svd) #$d is singular values, $u is vector of shape weights, $vt is vector of gene weights

gm_2bpls$svd$d[1]/sum(gm_2bpls$svd$d) #first singular value accounts for ... covariation

hist(abs(gm_2bpls$right.pls.vectors[,1]))#larger absolute value, the more weight they contribute
which.max(abs(gm_2bpls$right.pls.vectors[,1]))
gm_2bpls$right.pls.vectors[175,1]

sort(abs(gm_2bpls$right.pls.vectors[,1]), decreasing = TRUE)
summary(gm_2bpls)

#### sparse PLS regression with mixOmics ####
shape_subset_mat <- two.d.array(shape_subset_arr)
X <- shape_subset_mat
Y <- gene_subset

spls_results <- spls(X, Y, keepX=c(96,96), keepY = c(5,5)) #Keep all landmark coordinates but only top 5 genes
gene_loadings <- spls_results$loadings$Y #give the loading vectors their own object 
head(sort(abs(gene_loadings[,1]), decreasing = TRUE)) #sort the loading vectors of the first component 
head(sort(abs(gene_loadings[,2]), decreasing = TRUE)) #sort the loading vectors of the second component

#### PLS regression (SIMPLS) with pls ####
simpls_pls <- simpls.fit(X, Y, ncomp = 2) #fit the regression using the SIMPLS method
head(sort(abs(simpls_pls$Yloadings[,1]), decreasing = TRUE)) #sort the loading vectors of the first component 
head(sort(abs(simpls_pls$Yloadings[,2]), decreasing = TRUE)) #sort the loading vectors of the second component

