#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_filtRawBam
##SBATCH --output=ATACpipelineDev_filtRawBam.%a.out
##SBATCH --error=ATACpipelineDev_filtRawBam.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=8:00:00
##SBATCH --array=0-72
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA734635/*raw.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_bam} | sed 's/.raw.bam$//g'`

module load samtools

# keep proper pair
# filter out mate unmapped
# filter out poor quality
#-F 516 is the single end tag, 524 is the paired end tag

#samtools view -bh -F 516 -q 30 ${this_bam} > ${prefix}.filt.bam
samtools view -bh -f 2 -F 524 -q 30 SRR15054980.raw.bam > SRR15054980.filt.bam
#sort by read name, fixmate for marking duplicates
#samtools sort -n ${prefix}.filt.bam > ${prefix}.filt.qnameSort.bam
samtools sort -n SRR15054980.filt.bam > SRR15054980.filt.qnameSort.bam

#samtools fixmate ${prefix}.filt.qnameSort.bam ${prefix}.fixmate.bam
samtools fixmate SRR15054980.filt.qnameSort.bam SRR15054980.fixmate.bam

#1796 is the single end tag, 1804 is the paired end tag
#samtools view -bh -F 1796 ${prefix}.fixmate.bam > ${prefix}.fixmate.filt.bam 
samtools view -bh -f 2 -F 1804 SRR15054980.fixmate.bam > SRR15054980.fixmate.filt.bam 

#samtools sort ${prefix}.fixmate.filt.bam > ${prefix}.toMark.bam
samtools sort SRR15054980.fixmate.filt.bam > SRR15054980.toMark.bam

rm SRR15054980.fixmate.filt.bam

rm SRR15054980.fixmate.bam

rm SRR15054980.filt.qnameSort.bam

rm SRR15054980.filt.bam

