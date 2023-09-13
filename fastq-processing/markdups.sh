#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_markDup
#SBATCH --output=ATACpipelineDev_markDup.%a.out
#SBATCH --error=ATACpipelineDev_markDup.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1:00:00
#SBATCH --array=215-216
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*toMark.bam)

this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

prefix=`echo ${this_bam} | sed 's/.toMark.bam$//g'`

in_picard=/Genomics/grid/users/alea/programs/picard-tools-1.141/picard.jar

in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0

module load R
module load java
module load samtools
module load GATK/4.3.0.0

echo $prefix
echo $this_bam

gatk  MarkDuplicates \
-I ${this_bam} \
-O ${prefix}.markedDup.bam \
-M ${prefix}.markedDup.metrics.txt \
--ASSUME_SORTED true \
--VALIDATION_STRINGENCY LENIENT

echo "done"
