#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_getFRiP
##SBATCH --output=ATACpipelineDev_FRiPscore.%a.out
##SBATCH --error=ATACpipelineDev_FRiPscore.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
##SBATCH --array=0-215
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*no-blacklist.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_bam} | sed 's/.no-blacklist.bam$//g'`

#peakfiles=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/MACS2/SRR*_peaks.saf)
#this_peakfile=${peakfiles[${SLURM_ARRAY_TASK_ID}]}
#prefix2=`echo ${this_peakfile} | sed 's/_peaks.saf$//g'`


#Assign Reads to Peaks
#featureCounts -a ${this_peakfile} -F SAF -B -o ${prefix2}_FRiP.txt -p --largestOverlap -T 4 ${this_bam}

#Create File Listing FRiP for All Samples
touch FRiP_scores.txt

samples='/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/MACS2/ATACpipelineDev_FRiPscore.*.err'

for sample in ${samples}
do
grep "Load annotation file" ${sample} >> FRiP_scores.txt
done


touch temp.txt

for sample in ${samples}
do
grep "Successfully assigned alignments" ${sample} >> temp.txt
done

paste FRiP_scores.txt temp.txt > tmp
mv tmp FRiP_scores.txt
