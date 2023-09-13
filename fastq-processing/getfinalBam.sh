#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_finalBam
#SBATCH --output=ATACpipelineDev_finalBam.%a.out
#SBATCH --error=ATACpipelineDev_finalBam.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=215
#SBATCH --time=1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*dedup.bam)

this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

prefix=`echo ${this_bam} | sed 's/.dedup.bam$//g'`

module load samtools

CHROMOSOMES=$(samtools view -H ${this_bam} | grep '^@SQ' | cut -f2 | grep -v -e Y -e X -e M -e gl -e _random | sed 's/SN://' | xargs echo)

#CHROM=$(samtools view -H ${this_bam} | grep '^@SQ' | cut -f2 | grep -v -e Y -e X -e M -e gl -e _random | sed 's/SN:chr//' | xargs echo)

#echo $CHROMOSOMES
prefix=`echo ${this_bam} | sed 's/.filt.dedup.bam$//'`

## filter out all chrX/Y, mtDNA, and scaffold alignments

samtools view -bh -o ${prefix}.final.bam ${this_bam} $CHROMOSOMES

samtools index ${prefix}.final.bam

samtools idxstats ${prefix}.final.bam > ${prefix}.final.bam.idxstats
