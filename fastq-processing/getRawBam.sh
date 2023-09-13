#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_getBam
##SBATCH --output=ATACpipelineDev_getBam.%a.out
##SBATCH --error=ATACpipelineDev_getBam.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
##SBATCH --array=0-72
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

#sams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/SRR15054980.trimmed.sam)

#this_sam=${sams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_sam} | sed 's/.sam//g'`

module load samtools

samtools view -h -b /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/SRR15054980.trimmed.sam | samtools sort -o SRR15054980.raw.bam - 

samtools index SRR15054980.raw.bam

#get stats for mtDNA and other alignment info
samtools idxstats SRR15054980.raw.bam > SRR15054980.raw.bam.stats


