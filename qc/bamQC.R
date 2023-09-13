#!/bin/bash
#SBATCH --job-name=bamQC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

library(GenomicRanges)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ChIPpeakAnno)
library(GenomicAlignments)

bamfile <- "/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/SRR15054764.final.bam"
bamQC(bamfile, index = bamfile, outPath = "/Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/SRR15054764.final.bamqc.bam"))
