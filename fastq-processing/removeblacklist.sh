#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_removeblacklist
#SBATCH --output=ATACpipelineDev_removeblacklist.%a.out
#SBATCH --error=ATACpipelineDev_removeblacklist.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00
#SBATCH --array=215
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*final.bam)

this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

prefix=`echo ${this_bam} | sed 's/.final.bam$//g'`

#blacklist=/Genomics/ayroleslab2/emma/ATAC-Seq/hg19-blacklist.v2.bed.gz

module load samtools
module load bedtools 

#Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -a ${this_bam} -b /Genomics/ayroleslab2/emma/ATAC-Seq/hg19-blacklist.v2.bed.gz > ${prefix}.no-blacklist.bam

