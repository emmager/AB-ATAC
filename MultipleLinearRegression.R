library(ggplot2)
library(ggbiplot)
library(lme4)
library(variancePartition)
library(doParallel)
library(foreach)
library(patchwork)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(mixOmics)

setwd("/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248")
mydir <- "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/"
higherdir <- "/Genomics/ayroleslab2/emma/ATAC-Seq/"

load( "temp.RData")

source( paste0( scriptDir, "myFunction_fxn.R"))

dplyr::select()

#Check if reads from all bins in 1 ind follow a normal distribution
SRR15054764_binnedreads <- read.csv(file = paste0(mydir, "SRR15054764/SRR15054764.100kbCov_final.txt"), head = FALSE, sep = "\t")
SRR15054764_binnedreads$log2transformed_TPM <- log2(SRR15054764_binnedreads$V9 + .01)
within_sample_TPMs <- hist(SRR15054764_binnedreads$log2transformed_TPM) #long left tail

png(file= paste0(mydir,"within_sample_TPMs.png")
hist(SRR15054764_binnedreads$log2transformed_TPM, breaks=100)
dev.off()

#Check if the reads in 1 bin from all inds follow a normal distribution
all_bin_TPMs <- read.csv(file = paste0(mydir,"perbinTPMs.txt", head = TRUE, sep = "\t", row.names= 1))
all_bin_TPMs <- t(all_bin_TPMs)
all_bin_TPMs <- as.data.frame(all_bin_TPMs)
bin1_TPMs <- as.data.frame(all_bin_TPMs$Bin8) #rows = individuals, col= bin 1 TPM
colnames(bin1_TPMs) <- "Bin1"
bin1_TPMs$log2transformed <- log2(bin1_TPMs$Bin1 + .01)

within_bin_TPMs <- hist(bin1_TPMs$log2transformed)

png(file= paste0(mydir, "within_bin_TPMs.png")
hist(bin1_TPMs$log2transformed, breaks=100)
dev.off()


#Violin Plot to determine important factors

##Import Read Data
rowstoremove <- read.csv(file = paste0(higherdir,"removelines.txt"), head = FALSE, sep = "\t")
log2all_bin_TPMs <- read.csv(file = paste0(mydir, "perbinTPMs.txt"), head= TRUE, sep = "\t", row.names= 1)
log2all_bin_TPMsfilt <- log2all_bin_TPMs[-c(rowstoremove$V1),-171] #remove the last column which is all NAs and rows which are poorly annotated/telomeres.
#write.csv(log2all_bin_TPMsfilt, file= paste0(mydir, "log2all_bin_TPMsfilt.txt")

##Import Metadata
HQ_sample_metadata <- read.csv(file = paste0(mydir, "highqualsamples_metadata.txt"), head = TRUE, sep = "\t", row.names= 1)


#Formula indicating which Variables to include in analysis
formula <- ~ (1|Cell_Subset) + (1|Donor_ID) + (1|Sex) + (1|Age_Group) + FRiP_Score + ATAC_Preparation_Group + Percent_Dups

browser()

varPart<-fitExtractVarPartModel(log2all_bin_TPMsfilt[1:26574,],formula,HQ_sample_metadata)

#sortvariables(i.e.columns)bymedianfractionofvarianceexplained 
vp<-sortCols(varPart)
plotVarPart(vp)


png(file= paste0(mydir, "covarexplained.png"))
plotVarPart(vp)
dev.off()

#Make sub files that are unique to each Cell Type
(unique(HQ_sample_metadata$Cell_Subset))
#metadata
for (i in 1:length(unique(HQ_sample_metadata$Cell_Subset))){
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_metadata"), HQ_sample_metadata[HQ_sample_metadata$Cell_Subset == unique(HQ_sample_metadata$Cell_Subset)[i],])
}
#perbinlog2BPMs
for (i in 1:length(unique(HQ_sample_metadata$Cell_Subset))){
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_BPMs"), log2all_bin_TPMsfilt[,c(rownames(HQ_sample_metadata[HQ_sample_metadata$Cell_Subset == unique(HQ_sample_metadata$Cell_Subset)[i],]))])
}

#Per Cell line Violin Plots
BPM_files <- (ls(pattern = "BPMs"))
metadata_files <- (ls(pattern = "metadata"))

formula <- ~ (1|Sex) + (1|Age_Group) + FRiP_Score + ATAC_Preparation_Group + Percent_Dups

browser()

for (i in 1:length(unique(HQ_sample_metadata$Cell_Subset))){
  varPart<-fitExtractVarPartModel((get(BPM_files[i]))[1:26574,],formula,(get(metadata_files[i])))
  vp<-sortCols(varPart)
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_variancePar"), plotVarPart(vp, label.angle = 60, main = unique(HQ_sample_metadata$Cell_Subset)[i]))
}

mylist <- list()
for (i in 1:length(unique(HQ_sample_metadata$Cell_Subset))){
  p <- get(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_variancePar"))
  mylist[[i]] <- p
 } 

 png(file= "CellTypeSpec_VariancePar.png")
 patchwork::wrap_plots(mylist, nrow = 4, ncol = 4)
dev.off()


#Models
#Make table with columns = sample ID and rows being per bin counts and metadata
df <- HQ_sample_metadata %>% separate(Cell_Subset, c('CD', 'Subtype'))

#pca before 
forpca <- prcomp( t( log2all_bin_TPMsfilt ), center = TRUE, scale = TRUE)
summary( forpca)

for (i in 1:length(HQ_sample_metadata)){
png(paste0(mydir,colnames(HQ_sample_metadata[i]),"_prepca.png"))
print((ggbiplot( forpca, var.axes = F) + geom_point( aes( color = HQ_sample_metadata[,i]))))
dev.off()
}

#Determine which factors are contributing to the principal components
HQ_sample_metadata$Dummy_Sex <- ifelse(HQ_sample_metadata$Sex == "Female", 1, 0)
HQ_sample_metadata$Dummy_AgeGroup <- ifelse(HQ_sample_metadata$Age_Group == "Young", 1, 0)

techFactorMatrix <- cbind( HQ_sample_metadata[,c(4:6,8,10:11)], forpca$x[ , 1:10])
mycor <- cor( techFactorMatrix, method = "spearman")

png(file= paste0(mydir,"2023-09-16_TechFactorCorr_PCA.png"))
ggcorrplot(mycor)
dev.off()

fviz_contrib(forpca, choice = "var", axes = 1, top = 10)
#Trying to make column w/ CD4, CD8 or CD8 Naive
df <- HQ_sample_metadata %>% separate(Cell_Subset, c('CD', 'Subtype'))
HQ_sample_metadata$Immune_Cell_Type <- df$CD
(ggbiplot( forpca, var.axes = F) + geom_point( aes( color = df$CD)))






ggcorrplot( mycor)

corr_matrix <- cor(forpca)
ggcorrplot(corr_matrix)

#data.pca <- PCA(corr_matrix, graph = TRUE)
#summary(data.pca)

#plot(data.pca, choix = "ind", axes= c(1,2))

#fviz_screeplot(data.pca, ncp=10)

#registerDoParallel(cl)
residual <- matrix(0, nrow = 26574, ncol = 170)
for (i in 1:nrow(log2all_bin_TPMsfilt)){
  mod <- lm(unlist(log2all_bin_TPMsfilt[i,]) ~  Cell_Subset + Age_Group + FRiP_Score + Sex, data = HQ_sample_metadata)
  residual[i,] <- residuals(mod)
}

#Check if residuals from all bins in 1 ind follow a normal distribution
residual <- as.data.frame(residual)
within_sample_resids <- hist(residual$V1, breaks = 100)

png(file= paste0(mydir,"within_sample_resids.png"))
within_sample_resids
dev.off()

#Check if the reads in 1 bin from all inds follow a normal distribution
within_bin_resids <- hist(unlist(residual[5,]), breaks = 100) #nottttttt reallly. Look into this. 

png(file= paste0(mydir,"within_bin_resids.png"))
within_bin_resids
dev.off()


#pca after
postpca <- prcomp( t( residual ), center = TRUE, scale = TRUE)
summary( postpca)

for (i in 1:length(HQ_sample_metadata)){
png(paste0(mydir,colnames(HQ_sample_metadata[i]),"_postpca.png"))
print((ggbiplot( postpca, var.axes = F) + geom_point( aes( color = HQ_sample_metadata[,i]))))
dev.off()
}



#Correct all Cell Type Specific TPMs for technical factors

for (i in 1:length(unique(HQ_sample_metadata$Cell_Subset))){
  residual <- matrix(0, nrow = nrow((get(BPM_files[i]))), ncol = ncol((get(BPM_files[i]))))
  for (j in 1:nrow((get(BPM_files[i])))){
    mod <- lm(unlist(get(BPM_files[i])[j,]) ~ FRiP_Score + Sex, data = (get(metadata_files[i])))
    residual[j,] <- residuals(mod)
  }
  residual <- as.data.frame(residual)
  rownames(residual) <- rownames(log2all_bin_TPMsfilt)
  colnames (residual) <- colnames(get(BPM_files[i]))
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_residual"),residual)
}


residual_files <- (ls(pattern = "_residual"))
perchrbins <- read.csv(file = paste0(higherdir, "Chr_Bin_Boundaries.txt"), head = TRUE, sep = "\t") %>% unique
perchrbinsfilt <- perchrbins[-c(rowstoremove$V1),]


for (i in 1:12){
  pca <- prcomp( t(get(residual_files[i])), center = TRUE, scale = TRUE)
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[i],"_postpca"), pca)
}

pca_files <- (ls(pattern = "_postpca"))

for (j in 1:length(pca_files)){
  for (i in 1:length(get(metadata_files[j]))){
    png(paste0(mydir, unique(HQ_sample_metadata$Cell_Subset)[j] , "_", colnames(HQ_sample_metadata[i]), "_postpca.png"))
    print(ggbiplot((get(pca_files[j])), var.axes = F) + geom_point( aes( color = (get(metadata_files[j])[,i]))))
    dev.off()
}
}

#Make Correlation Matrix of Residuals

#Pull out the chr specific rows
for (i in 1:22){
for (j in 1:length(residual_files)){
  chr_resid <- get(residual_files[j])[c(perchrbins[perchrbins$Chr== paste0("chr",i),1]),]
  chr_resid <- na.omit(chr_resid)
  assign(paste0("chr",i,"_",unique(HQ_sample_metadata$Cell_Subset)[j],"_residual"), chr_resid)
}
}

chr_resid_files <- (ls(pattern = "chr"))
chr_resid_files <- (unique(chr_resid_files[ which( !( grepl( "bins", chr_resid_files)))]))[-1]

#Check for Nas in chr_resid_files
stand <- matrix(0, nrow = 2224, ncol = 1)
for (i in 1:nrow(get(chr_resid_files[2]))){
  stand[i,] <- sd(get(chr_resid_files[2])[i,])
}


cormatr <-  function(x) {
    return(cor((t(get(x)), method = "pearson") ))
}

#for (i in length(chr_resid_files)) {
cell_cors <- cor((t(get(chr_resid_files[1]))), method = "pearson")
assign(paste0(chr_resid_files[1],"_cor"), cell_cors)
}



#pearson vs. spearman?

cor_plots <- ls(pattern = "_cor")
cor_plots <- cor_plots[ which( ( grepl( "residual", cor_plots)))]

#png(paste0(mydir,"chr22_CD4_BulkNonNaive_corplots.png"))
ggcorrplot(chr1_CD4_BulkNaive_residual_cor)
#dev.off()


#Call First Eigenvetor with mixOmics

for (i in 1:length(cor_plots)){
  eigen <- nipals(get(cor_plots[i]), ncomp = 10)
 assign(paste0(cor_plots[i],"_eigenvector"), eigen)
}

eigen <- nipals(get(cor_plots[1]), ncomp = 1)
chrBinEigens <- data.frame(binCorrEigen = eigen$t)
head(chrBinEigens)
chrBinEigens$PC1

#Visualize Eigenvectors
# Read the data from the text file:
abRead <- function(file, chr=NULL){
	if (is.null(chr)){
		data <- read.table(file, head=TRUE)
	} else {
		data <- read.table(file, head=TRUE)
		data <- data[data$chr==chr,]
	}
	data
}

# barplot with nice colors:
abBarplot <- function(abdata, main = "", ylim = c(-2, 2), 
    top.col = "firebrick", bot.col = "grey50") {
	x <- abdata$eigen
    x <- as.numeric(x)
    n <- length(x)
    col <- rep(top.col, n)
    col[x < 0] <- bot.col
    barplot(x, ylim = ylim, bty = "n", xlab = "", ylab = "", 
        border = col, col = col, main = main, yaxt = "n")
}

abBarplot(chrBinEigens[1])

# convert the data to a GRanges:
convert2GRanges <- function(abdata){
	library(GenomicRanges)
	gr <- GRanges(seqnames=abdata$chr,IRanges(start=data$start, end=data$end))
	genome(gr) <- "hg19"
	gr$eigen <- abdata$eigen
	gr$domain <- abdata$domain
	gr
}


file <- "../data/prad_tumor_compartments_100kb.txt" # Should be path to the file
data <- abRead(file, "chr14")
abBarplot(data)

# To convert to a GRanges:
gr <- convert2GRanges(data)



















# Make functions
pcas <- function(x){
  pca <- prcomp( t(x)), center = TRUE, scale = TRUE)
  assign(paste0(unique(HQ_sample_metadata$Cell_Subset)[x],"_postpca"), pca)
}

#shapiro test and p value of esitmated effect size, Intagrative genomics. 


save.image( "temp.RData")


cell_cors <- cor((t(get(chr_resid_files[12]))), method = "pearson")
