#!/bin/bash

#SBATCH --job-name=num_of_reads
#SBATCH --output=PRJNA744248_allsamples_readcounts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --mail-type=end

module load samtools 

for file in /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*final.bam
do
prefix=`echo $file | sed 's/.final.bam$//'`
echo ${prefix}
samtools view  $file | wc -l 
done
