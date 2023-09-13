#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_dedup
#SBATCH --output=ATACpipelineDev_dedup.%a.out
#SBATCH --error=ATACpipelineDev_dedup.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00
#SBATCH --array=215
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*markedDup.bam)

this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

prefix=`echo ${this_bam} | sed 's/.markedDup.bam$//g'`

module load samtools


# keep proper pair
# filter out duplicates
# filter out mate unmapped
# filter out non-primary alignments
#Add -f 2 for paired end reads and adjust -F tag as noted in filterRawbam.sh

samtools view -bh -f 2 -F 1804 ${this_bam} > ${prefix}.filt.dedup.bam

samtools index ${prefix}.filt.dedup.bam

# final clean alignment stats

samtools idxstats ${prefix}.filt.dedup.bam > ${prefix}.stats
