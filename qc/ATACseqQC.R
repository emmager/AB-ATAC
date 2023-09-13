#Load Prerequisits
library(GenomicRanges)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ChIPpeakAnno)
library(GenomicAlignments)
library(phastCons100way.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(GenomicRanges)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ChIPpeakAnno)
library(GenomicAlignments)
library(phastCons100way.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(foreach)
library(doParallel)

print("all done loading packages")
args = commandArgs(trailingOnly=TRUE)
num <- as.numeric(args[1])

markedDup_Files <- list.files(path = "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248", pattern = "*.markedDup.bam$", full.names = TRUE)
final_Files <- list.files(path = "/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248", pattern = "*.no-blacklist.bam$", full.names = TRUE)

unfiltered_bamFileLabels <- gsub(".markedDup.bam","", basename(markedDup_Files[num]))
final_bamFileLabels <- gsub(".no-blacklist.bam","", basename(final_Files[num]))

#Library Complexity
## %MT reads, Dup rate, NRF, PBC_1, PBC_2, idxstats

this_unfiltered_bamQC <- bamQC(markedDup_Files[num])
unfiltered_stats <- this_unfiltered_bamQC[1:10]
write.csv(unfiltered_stats,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",unfiltered_bamFileLabels,"_unfiltered","_bamQC.txt"))

this_filtered_bamQC <- bamQC(final_Files[num])
filtered_stats <- this_filtered_bamQC[1:10]
write.csv(filtered_stats,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_final","_bamQC.txt"))

##Fragment Size Distribution

png(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",unfiltered_bamFileLabels,"_unfiltered","_fragSize.png"))
this_unfiltered_fragSize <- fragSizeDist(markedDup_Files[num],unfiltered_bamFileLabels)
dev.off()
write.csv(this_unfiltered_fragSize,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",unfiltered_bamFileLabels,"_unfiltered","_fragSize"))


png(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_final","_fragSize.png"))
this_filtered_fragSize <- fragSizeDist(final_Files[num],final_bamFileLabels)
dev.off()
write.csv(this_unfiltered_fragSize,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC",final_bamFileLabels,"_final","_fragSize"))


##Estimate Library Complexity

png(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",unfiltered_bamFileLabels,"_unfiltered","_libcomp.png"))
this_unfiltered_LibComplexity <- estimateLibComplexity(readsDupFreq(markedDup_Files[num]))
dev.off()
write.csv(this_unfiltered_LibComplexity,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",unfiltered_bamFileLabels,"_unfiltered","_libcomp"))

png(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_final","_libcomp.png"))
this_filtered_LibComplexity <- estimateLibComplexity(readsDupFreq(final_Files[num]))
dev.off()
write.csv(this_filtered_LibComplexity,file = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_filtered","_libcomp"))

#Peak Calling

## Shifting Aligned Reads

this_gal <- readBamFile(bamFile="/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/SRR15054980.no-blacklist.bam",tag=character(0), asMates=TRUE)
this_gal1 <- shiftGAlignmentsList(this_gal, outbam=paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_shifted.bam"))

#write.csv(this_gal1, file= paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_final","_shifted.bam"))

## Split Bam Files based on Nucleosome Free, Mono- di- or tri-nucleosomes
txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

obj_split <- splitGAlignmentsByCut(this_gal1,txs=txs, genome = "Hsapiens", outPath = paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_splitbam"))

dir(paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_splitbam")

#TSSE
this_shifted_bamfile <- paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_shifted.bam")
this_gal2<-readBamFile(bamFile=this_shifted_bamfile,tag=character(0), asMates=TRUE)
txs<-transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tsse<-TSSEscore(this_gal1,txs)
write.csv(tsse, file= paste0("/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/",final_bamFileLabels,"_final_tsse.txt")


#Plot Correlation??!


