#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_getbincov
##SBATCH --output=ATACpipelineDev_bincov.%a.out
##SBATCH --error=ATACpipelineDev_bincov.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
##SBATCH --array=0-216%80
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

##FeatureCounts Info
#	-p species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
#	-O assigns reads to all their overlapping meta-features.
#	-T specifies the number (n) of threads to be used.
#	-a is the genome annotation file (example_genome_annotation.gtf).
#	-o specifies the name of the output file, which includes the read counts (example_featureCounts_output.txt).
#	sorted_example_alignment.bam is an alignment file: in this file, the reads we want to count are aligned to the same genome as the annotation file.

##Remember to activate conda environment "subread" in argo-comp1
#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*no-blacklist.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_bam} | sed 's/.no-blacklist.bam$//g'`

featureCounts -a /Genomics/ayroleslab2/emma/ATAC-Seq/hg19_100kbBins_noBlacklist.saf -F SAF -B -o SRR15054980_100kbBinCoverage.txt -p --largestOverlap -T 4 SRR15054980.no-blacklist.bam
