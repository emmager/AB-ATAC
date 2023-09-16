library(ggplot2)
library(ggbiplot)
library(lme4)
library(tidyr)
library(dplyr)
library(variancePartition)
library(doParallel)
library(foreach)
library(patchwork)
library(corrr)
library(ggcorrplot)
library(FactoMineR)


setwd("/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248")


#Check if reads from all bins in 1 ind follow a normal distribution
SRR15054764_binnedreads <- read.csv(file ="/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/SRR15054764/SRR15054764.100kbCov_final.txt", head = FALSE, sep = "\t")
SRR15054764_binnedreads$log2transformed_TPM <- log2(SRR15054764_binnedreads$V9 + .01)
within_sample_TPMs <- hist(SRR15054764_binnedreads$log2transformed_TPM) #long left tail

png(file= "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/within_sample_TPMs.png")
hist(SRR15054764_binnedreads$log2transformed_TPM, breaks=100)
dev.off()

#Check if the reads in 1 bin from all inds follow a normal distribution
all_bin_TPMs <- read.csv(file ="/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/perbinTPMs.txt", head = TRUE, sep = "\t", row.names= 1)
all_bin_TPMs <- t(all_bin_TPMs)
all_bin_TPMs <- as.data.frame(all_bin_TPMs)
bin1_TPMs <- as.data.frame(all_bin_TPMs$Bin8) #rows = individuals, col= bin 1 TPM
colnames(bin1_TPMs) <- "Bin1"
bin1_TPMs$log2transformed <- log2(bin1_TPMs$Bin1 + .01)

within_bin_TPMs <- hist(bin1_TPMs$log2transformed)

png(file= "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/within_bin_TPMs.png")
hist(bin1_TPMs$log2transformed, breaks=100)
dev.off()


#Violin Plot to determine important factors

##Import Read Data
rowstoremove <- read.csv(file = "/Genomics/ayroleslab2/emma/ATAC-Seq/removelines.txt", head = FALSE, sep = "\t")
log2all_bin_TPMs <- read.csv(file ="/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/perbinTPMs.txt", head= TRUE, sep = "\t", row.names= 1)
log2all_bin_TPMsfilt <- log2all_bin_TPMs[-c(rowstoremove$V1),-171] #remove the last column which is all NAs and rows which are poorly annotated/telomeres.https://vscode-remote+ssh-002dremote-002bargo-002dcomp1-002eprinceton-002eedu.vscode-resource.vscode-cdn.net/tmp/RtmpNi8grw/vscode-R/plot.png?version%3D1694830401398

##Import Metadata
HQ_sample_metadata <- read.csv(file ="/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/highqualsamples_metadata.txt", head = TRUE, sep = "\t", row.names= 1)


#Formula indicating which Variables to include in analysis
formula <- ~ (1|Cell_Subset) + (1|Donor_ID) + (1|Sex) + (1|Age_Group) + FRiP_Score + ATAC_Preparation_Group + Percent_Dups

browser()

varPart<-fitExtractVarPartModel(log2all_bin_TPMsfilt[1:26574,],formula,HQ_sample_metadata)

#sortvariables(i.e.columns)bymedianfractionofvarianceexplained 
vp<-sortCols(varPart)
plotVarPart(vp)


png(file= "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/covarexplained.png")
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
png(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/",colnames(HQ_sample_metadata[i]),"_prepca.png"))
print((ggbiplot( forpca, var.axes = F) + geom_point( aes( color = HQ_sample_metadata[,i]))))
dev.off()
}
techFactorMatrix <- cbind( techFactorMatrix, forpca$x[ , 1:10])

#Trying to make column w/ CD4, CD8 or CD8 Naive
df <- HQ_sample_metadata %>% separate(Cell_Subset, c('CD', 'Subtype'))
HQ_sample_metadata$Immune_Cell_Type <- df$CD
(ggbiplot( forpca, var.axes = F) + geom_point( aes( color = df$CD)))




mycor <- cor( techFactorMatrix, method = "spearman")

ggcorrplot( mycor)

corr_matrix <- cor(forpca)
ggcorrplot(corr_matrix)

#data.pca <- PCA(corr_matrix, graph = TRUE)
#summary(data.pca)

#plot(data.pca, choix = "ind", axes= c(1,2))

#fviz_screeplot(data.pca, ncp=10)

#registerDoParallel(cl)
residual <- matrix(0, nrow = 26574, ncol = 170)
# for (i in 1:nrow(log2all_bin_TPMsfilt)){
  for (i in 1:10){
  mod <- lm(unlist(log2all_bin_TPMsfilt[i,]) ~  Cell_Subset + Age_Group + FRiP_Score + Sex, data = HQ_sample_metadata)
  residual[i,] <- residuals(mod)
}


residual <- matrix(0, nrow = 26574, ncol = 170)
# for (i in 1:nrow(log2all_bin_TPMsfilt)){
  for (i in 1:10){
  mod <- lm(unlist(log2all_bin_TPMsfilt[i,]) ~  Cell_Subset + Age_Group + FRiP_Score + Sex, data = HQ_sample_metadata)
  residual[i,] <- residuals(mod)
}

#Check if residuals from all bins in 1 ind follow a normal distribution
residual <- as.data.frame(residual)
within_sample_resids <- hist(residual$V1, breaks = 100)

png(file= "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/within_sample_resids.png")
within_sample_resids
dev.off()

#Check if the reads in 1 bin from all inds follow a normal distribution
within_bin_resids <- hist(unlist(residual[5,]), breaks = 100) #nottttttt reallly. Look into this. 

png(file= "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/within_bin_resids.png")
within_bin_resids
dev.off()

#shapiro test and p value of esitmated effect size, Intagrative genomics. 